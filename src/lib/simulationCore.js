// simulationCore.js
// ------------------------------------------------------------
// Mechanik (Rad):
//   I_disc * dω/dt = τ_drive - τ_brake
//   dθ/dt = ω
//
// Szenarien:
//   - "Constant Speed":
//       Motor kompensiert τ_brake → ω ~ konstant, Bremse heizt dauerhaft.
//   - "Emergency Brake":
//       τ_drive = 0, Abbremsen nur über Bremse (inkl. Haftreibung).
//   - "Acceleration with Brake":
//       τ_drive = τ_brake + τ_extra  (τ_extra aus fester A_POS_AWB)
//       → Rad beschleunigt trotz Bremse, Bremse heizt dauerhaft.
//
// Thermik:
//   - pulseWheelSeq[i] = Energie, die im Intervall i*dt..(i+1)*dt
//     in den betrachteten Radpunkt geht [J]
//   - 3D-Green-Funktion (Halbraum, surfaceFactor=2), Singularitätsfix
//     r_eff^2 = dist^2 + sourceRadius^2
//   - 2-Körper-Split via Effusivität oder heatSplitWheel
//   - kLoss als exp(-lambdaLoss * tau) im Kernel
// ------------------------------------------------------------

function fft_inplace(re, im, invert = false) {
  const n = re.length;
  let j = 0;
  for (let i = 1; i < n; i++) {
    let bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }

  for (let len = 2; len <= n; len <<= 1) {
    const ang = (2 * Math.PI / len) * (invert ? -1 : 1);
    const wlenRe = Math.cos(ang);
    const wlenIm = Math.sin(ang);

    for (let i = 0; i < n; i += len) {
      let wRe = 1;
      let wIm = 0;

      for (let k = 0; k < len / 2; k++) {
        const uRe = re[i + k];
        const uIm = im[i + k];

        const vRe = re[i + k + len / 2] * wRe - im[i + k + len / 2] * wIm;
        const vIm = re[i + k + len / 2] * wIm + im[i + k + len / 2] * wRe;

        re[i + k] = uRe + vRe;
        im[i + k] = uIm + vIm;
        re[i + k + len / 2] = uRe - vRe;
        im[i + k + len / 2] = uIm - vIm;

        const nextWRe = wRe * wlenRe - wIm * wlenIm;
        const nextWIm = wRe * wlenIm + wIm * wlenRe;
        wRe = nextWRe;
        wIm = nextWIm;
      }
    }
  }

  if (invert) {
    for (let i = 0; i < n; i++) {
      re[i] /= n;
      im[i] /= n;
    }
  }
}

function convolveFFT_real(a, b) {
  const outLen = a.length + b.length - 1;
  let n = 1;
  while (n < outLen) n <<= 1;

  const are = new Float64Array(n);
  const aim = new Float64Array(n);
  const bre = new Float64Array(n);
  const bim = new Float64Array(n);

  are.set(a);
  bre.set(b);

  fft_inplace(are, aim, false);
  fft_inplace(bre, bim, false);

  const cre = new Float64Array(n);
  const cim = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    cre[i] = are[i] * bre[i] - aim[i] * bim[i];
    cim[i] = are[i] * bim[i] + aim[i] * bre[i];
  }

  fft_inplace(cre, cim, true);
  return cre.subarray(0, outLen);
}

function wrapAngle(x) {
  const TWO_PI = 2 * Math.PI;
  if (!Number.isFinite(x)) return 0;
  x = x % TWO_PI;
  if (x > Math.PI) x -= TWO_PI;
  if (x < -Math.PI) x += TWO_PI;
  return x;
}

function effusivity(k, rho, cp) {
  const kk = Math.max(1e-12, k);
  const rr = Math.max(1e-9, rho);
  const cc = Math.max(1e-9, cp);
  return Math.sqrt(kk * rr * cc);
}

export function simulationCore(params) {
  // dt ist rein numerisch, nicht im UI
  const dt = 1e-3; // 1 ms
  const duration = Math.max(dt, params.duration || 5);
  const steps = Math.max(1, Math.floor(duration / dt) + 1);

  // Feste interne Beschleunigung für AWB / Emergency
  const A_POS_AWB = 0.5;     // m/s²  (Anfahren gegen Bremse)
  const A_NEG_EMERGENCY = -1.5; // m/s² (für Notbremsen falls Motorbremse dazu käme)

  // --- Material Rad ---
  const rho_w = Number.isFinite(params.rho_disc) ? Math.max(1e-9, params.rho_disc) : 7800;
  const cp_w  = Number.isFinite(params.cp_disc)  ? Math.max(1e-9, params.cp_disc)  : 460;
  const k_w   = Number.isFinite(params.k_disc)   ? Math.max(1e-12, params.k_disc)  : 50;
  const alpha_w = k_w / (rho_w * cp_w);

  // --- Material Schiene / Umgebung (für Split) ---
  const rho_r = Number.isFinite(params.rho_rail) ? Math.max(1e-9, params.rho_rail) : rho_w;
  const cp_r  = Number.isFinite(params.cp_rail)  ? Math.max(1e-9, params.cp_rail)  : cp_w;
  const k_r   = Number.isFinite(params.k_rail)   ? Math.max(1e-12, params.k_rail)  : k_w;

  const T0_C = Number.isFinite(params.T0) ? params.T0 : 20;

  // Messpunkt-Abstand (oder Tiefe)
  const dist = Number.isFinite(params.distance) ? Math.max(0, params.distance) : 0.001;
  const distSq = dist * dist;

  // Szenario
  const scenario = typeof params.scenario === "string" ? params.scenario : "Constant Speed";
  const phi_obs = Number.isFinite(params.phi_obs) ? params.phi_obs : 0;

  // Rad-Geometrie
  const R = Number.isFinite(params.R) ? Math.max(1e-9, params.R) : 0.3;

  // Bremsklotz-Geometrie
  const padLength = Number.isFinite(params.padLength) ? Math.max(1e-6, params.padLength) : 0.05;
  const maxPadLength = Math.PI * R;
  const padLengthClamped = Math.min(padLength, maxPadLength);

  const contactAngle = padLengthClamped / R;
  const contactArc   = Math.max(1e-12, padLengthClamped);

  const padWidth = Number.isFinite(params.padWidth) ? Math.max(1e-6, params.padWidth) : 0.05;
  const contactArea = padLengthClamped * padWidth;

  const contactDepth = Number.isFinite(params.contactDepth)
    ? Math.max(1e-6, params.contactDepth)
    : 0.002;

  const thermalVolumeFactor = Number.isFinite(params.thermalVolumeFactor)
    ? Math.max(1e-6, params.thermalVolumeFactor)
    : 1.0;

  const V_effect = contactArea * contactDepth * thermalVolumeFactor;

  // Trägheitsmoment Rad
  let I_disc = Number.isFinite(params.I_disc) ? Math.max(1e-12, params.I_disc) : NaN;
  if (!Number.isFinite(I_disc)) {
    const wheelMass = Number.isFinite(params.wheelMass) ? Math.max(1e-9, params.wheelMass) : 200;
    I_disc = 0.5 * wheelMass * R * R;
  }

  // Reibung
  const mu = Number.isFinite(params.mu) ? Math.max(0, params.mu) : 0.4;

  // FN ist GESAMTNORMALKRAFT (für dieses Rad), nicht pro Pad
  const nBK_raw = Number.isFinite(params.nBK) ? params.nBK : 1;
  const nBK = Math.max(1, Math.min(2, Math.floor(nBK_raw)));
  const FN_total = Number.isFinite(params.FN) ? Math.max(0, params.FN) : 3000;
  const FN_pad = nBK > 0 ? FN_total / nBK : 0;   // Aufteilung auf Pads

  const mu_s = Number.isFinite(params.mu_s) ? Math.max(mu, params.mu_s) : mu;

  const frictionToHeat = Number.isFinite(params.frictionToHeat)
    ? Math.max(0, Math.min(1, params.frictionToHeat))
    : 0.9;

  // Singularitätsfix: effektiver Quellradius
  const baseRadius = padLengthClamped / 40;
  const sourceRadius = Number.isFinite(params.sourceRadius)
    ? Math.max(1e-4, params.sourceRadius)
    : Math.max(1e-4, baseRadius);

  const rEffSq = distSq + sourceRadius * sourceRadius;

  // 2-Körper-Split Rad/Schiene
  let heatSplitWheel = Number.isFinite(params.heatSplitWheel)
    ? Math.max(0, Math.min(1, params.heatSplitWheel))
    : NaN;

  if (!Number.isFinite(heatSplitWheel)) {
    const e_w = effusivity(k_w, rho_w, cp_w);
    const e_r = effusivity(k_r, rho_r, cp_r);
    heatSplitWheel = e_w / (e_w + e_r);
    if (!Number.isFinite(heatSplitWheel)) heatSplitWheel = 0.5;
  }
  const heatSplitRail = 1 - heatSplitWheel;

  // Wärmeabgabe / Konvektion (als Kernel-Dämpfung)
  const h_loss = Number.isFinite(params.kLoss) ? Math.max(0, params.kLoss) : 200;
  const Cw_local = Math.max(1e-12, rho_w * cp_w * V_effect);
  const lambdaLoss = (h_loss * contactArea) / Cw_local;

  // Halbraum-Faktor (Oberfläche)
  const surfaceFactor = 2;

  // Bremsklotz-Winkel
  let bkAngles = Array.isArray(params.bkAngles) && params.bkAngles.length > 0
    ? params.bkAngles
    : (nBK === 1 ? [0] : [0, Math.PI]); // z.B. gegenüberliegend

  // Startbedingungen ω
  const hasV0 = Number.isFinite(params.v0);
  const v0 = hasV0 ? Math.max(0, params.v0) : null;
  const wheelFrequency = Number.isFinite(params.wheelFrequency)
    ? Math.max(0, params.wheelFrequency)
    : 12;

  let omega;
  if (scenario === "Acceleration with Brake" && !hasV0) {
    omega = 0; // echtes Anfahren
  } else {
    omega = hasV0 ? v0 / R : 2 * Math.PI * wheelFrequency;
  }
  let v_mps = omega * R;

  const pulseWheelSeq = new Float64Array(steps);
  const pulseRailSeq  = new Float64Array(steps);
  const velocityProfile = new Float64Array(steps);

  let theta = 0;
  const omegaEps = 1e-6;
  const NsubBase = 8;

  const tauFromAccel = (aCmd) => I_disc * (aCmd / R);

  for (let i = 0; i < steps; i++) {
    // Substepping, um den Kontakt besser aufzulösen
    const dThetaStepEst = Math.abs(omega * dt);
    const Nsub = Math.max(
      NsubBase,
      Math.ceil(dThetaStepEst / Math.max(1e-12, contactAngle * 0.25))
    );
    const dtSub = dt / Nsub;

    for (let s = 0; s < Nsub; s++) {
      let tauBrakeKinetic = 0;
      let tauBrakeStaticMax = 0;
      let padsInContact = 0;

      // Beobachtungspunkt auf dem Rad
      const obsAngle = wrapAngle(theta + phi_obs);

      for (let j = 0; j < nBK; j++) {
        const bAngle = wrapAngle(
          Number.isFinite(bkAngles[j]) ? bkAngles[j] : 0
        );
        const d = Math.abs(wrapAngle(obsAngle - bAngle));
        const inContact = d <= contactAngle * 0.5;
        if (!inContact) continue;

        padsInContact++;

        // WICHTIG: FN_pad (Gesamt FN wird auf Pads verteilt)
        const tauPadKin = mu * FN_pad * R;
        const tauPadStaticMax = mu_s * FN_pad * R;

        tauBrakeKinetic += tauPadKin;
        tauBrakeStaticMax += tauPadStaticMax;
      }

      const inAnyContact = padsInContact > 0;

      // τ_drive je Szenario
      let tauDrive = 0;
      if (scenario === "Constant Speed") {
        // Motor kompensiert die Bremse → ω konstant
        tauDrive = inAnyContact ? tauBrakeKinetic : 0;
      } else if (scenario === "Emergency Brake") {
        // Kein Motorantrieb
        tauDrive = 0;
      } else if (scenario === "Acceleration with Brake") {
        // Motor: Bremse + zusätzliche Beschleunigung
        const tauExtra = tauFromAccel(A_POS_AWB);
        tauDrive = inAnyContact ? (tauBrakeKinetic + tauExtra) : tauExtra;
      }

      // Haftreibung nur für Emergency Brake:
      if (
        scenario === "Emergency Brake" &&
        inAnyContact &&
        Math.abs(omega) < omegaEps
      ) {
        const tauDrivePos = Math.max(0, tauDrive);
        if (tauDrivePos <= tauBrakeStaticMax) {
          omega = 0;
          v_mps = 0;
          // Kein Schlupf → keine Heizleistung
          continue;
        }
      }

      const tauBrake = inAnyContact ? tauBrakeKinetic : 0;
      const tauNet = tauDrive - tauBrake;
      const alphaRot = tauNet / I_disc;

      // einfache Zeitintegration
      omega = omega + alphaRot * dtSub;
      if (omega < 0) omega = 0;
      v_mps = omega * R;

      theta = theta + omega * dtSub;

      // Reibungswärme nur bei Schlupf
      if (inAnyContact && omega > omegaEps) {
        const P_fric = tauBrake * omega; // W
        const dE_total = frictionToHeat * P_fric * dtSub;

      // Kein zusätzlicher Slip-Faktor → dein Punkt repräsentiert den "Hotspot"
        const dE_local = dE_total * heatSplitWheel;

        pulseWheelSeq[i] += dE_local;

        pulseRailSeq[i]  += heatSplitRail  * dE_local;
      }
    }

    velocityProfile[i] = v_mps;
  }

  // Leistung fürs Plotten
  const powerSeq = new Float64Array(steps);
  for (let i = 0; i < steps; i++) {
    powerSeq[i] = pulseWheelSeq[i] / dt;
  }

  // Green-Kernel (Antwort auf 1 Joule im Radpunkt)
  const kernel = new Float64Array(steps);
  for (let i = 0; i < steps; i++) {
    const tau = (i + 0.5) * dt;
    const denom = Math.pow(4 * Math.PI * alpha_w * tau, 1.5);
    let G = 0;
    if (denom > 0) {
      G = Math.exp(-rEffSq / (4 * alpha_w * tau)) / (denom * rho_w * cp_w);
    }
    kernel[i] = surfaceFactor * G * Math.exp(-lambdaLoss * tau);
  }

  const conv = convolveFFT_real(pulseWheelSeq, kernel);

  // Ausgabe resamplen
  const maxPoints = Math.max(
    2,
    Math.min(5000, Number.isFinite(params.maxPoints) ? params.maxPoints : 500)
  );
  const Nout = Math.min(maxPoints, steps);

  const sampleAt = (arr, tIndexFloat) => {
    if (tIndexFloat <= 0) return arr[0];
    const i0 = Math.floor(tIndexFloat);
    const i1 = Math.min(arr.length - 1, i0 + 1);
    const w = tIndexFloat - i0;
    return (1 - w) * arr[i0] + w * arr[i1];
  };

  const out = [];
  for (let i = 0; i < Nout; i++) {
    const tFrac = Nout === 1 ? 0 : i / (Nout - 1);
    const tIndexFloat = tFrac * (steps - 1);
    const time = tFrac * duration;

    const pointTemp = T0_C + sampleAt(conv, tIndexFloat);

    // kumulative Energie im Rad
    const cumulativeEnergyWheel = (() => {
      const idx = Math.floor(tIndexFloat);
      let sum = 0;
      for (let j = 0; j <= idx; j++) sum += pulseWheelSeq[j];
      const frac = tIndexFloat - idx;
      if (idx + 1 < pulseWheelSeq.length) sum += pulseWheelSeq[idx + 1] * frac;
      return sum;
    })();

    out.push({
      time,
      pointTemp,
      cumulativeEnergy: cumulativeEnergyWheel,
      power: sampleAt(powerSeq, tIndexFloat),
      velocity: sampleAt(velocityProfile, tIndexFloat),
      phi_obs
    });
  }

  return out;
}
