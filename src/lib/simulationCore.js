// simulationCore.js
// ------------------------------------------------------------
// Dynamik:
//   I_eff * dω/dt = τ_drive - τ_brake
//   dθ/dt = ω
//
// I_eff = I_disc + mEff * R^2,  mEff = trainMass / nDrivenWheels (intern)
//
// Numerik:
// - Substepping + exaktes Update im Substep.
//
// Thermik:
// - pulseWheelSeq: Energie pro dt [J] (Rad-Anteil)
// - 3D-Green im Halbraum (surfaceFactor=2), Singularitäts-Fix r^2 + a^2
// - 2-Körper-Split via Effusivität oder heatSplitWheel
// - kLoss als exp(-lambda * tau) im Kernel
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
  const dt = Math.max(1e-6, params.dt || 1e-3);
  const duration = Math.max(dt, params.duration || 5);
  const steps = Math.max(1, Math.floor(duration / dt) + 1);

  // --- Material Rad ---
  const rho_w = Number.isFinite(params.rho_disc) ? Math.max(1e-9, params.rho_disc) : 7800;
  const cp_w  = Number.isFinite(params.cp_disc)  ? Math.max(1e-9, params.cp_disc)  : 460;
  const k_w   = Number.isFinite(params.k_disc)   ? Math.max(1e-12, params.k_disc)  : 50;
  const alpha_w = k_w / (rho_w * cp_w);

  // --- Material Schiene ---
  const rho_r = Number.isFinite(params.rho_rail) ? Math.max(1e-9, params.rho_rail) : rho_w;
  const cp_r  = Number.isFinite(params.cp_rail)  ? Math.max(1e-9, params.cp_rail)  : cp_w;
  const k_r   = Number.isFinite(params.k_rail)   ? Math.max(1e-12, params.k_rail)  : k_w;

  const T0_C = Number.isFinite(params.T0) ? params.T0 : 20;

  // Messabstand
  const dist = Number.isFinite(params.distance) ? Math.max(0, params.distance) : 0.001;
  const distSq = dist * dist;

  // Szenario
  const scenario = typeof params.scenario === "string" ? params.scenario : "Constant Speed";
  const aNeg = Number.isFinite(params.aNeg) ? params.aNeg : -0.8;
  const aPos = Number.isFinite(params.aPos) ? params.aPos : 0.4;
  const phi_obs = Number.isFinite(params.phi_obs) ? params.phi_obs : 0;

  // Geometrie
  const R = Number.isFinite(params.R) ? Math.max(1e-9, params.R) : 0.3;

  // Bremsklotz-Geometrie
  const padLength = Number.isFinite(params.padLength) ? Math.max(1e-6, params.padLength) : 0.05; // 50 mm
  const maxPadLength = Math.PI * R;
  const padLengthClamped = Math.min(padLength, maxPadLength);

  const contactAngle = padLengthClamped / R;
  const contactArc   = Math.max(1e-12, padLengthClamped);

  const padWidth = Number.isFinite(params.padWidth) ? Math.max(1e-6, params.padWidth) : 1.0;
  const contactArea = padLengthClamped * padWidth;

  const contactDepth = Number.isFinite(params.contactDepth)
    ? Math.max(1e-6, params.contactDepth)
    : 0.001;

  const thermalVolumeFactor = Number.isFinite(params.thermalVolumeFactor)
    ? Math.max(1e-6, params.thermalVolumeFactor)
    : 1.0;

  const V_effect = contactArea * contactDepth * thermalVolumeFactor;

  // Rad-Trägheitsmoment
  let I_disc = Number.isFinite(params.I_disc) ? Math.max(1e-12, params.I_disc) : NaN;
  if (!Number.isFinite(I_disc)) {
    const wheelMass = Number.isFinite(params.wheelMass) ? Math.max(1e-9, params.wheelMass) : 300;
    I_disc = 0.5 * wheelMass * R * R;
  }

  // --- Zug / Adhäsion (intern, keine Slider nötig) ---
  const g = 9.81;
  const mTrain = Number.isFinite(params.trainMass) ? Math.max(1, params.trainMass) : 40000; // kg
  const muAdh  = Number.isFinite(params.muAdh) ? Math.max(0, params.muAdh) : 0.15;
  const nDrivenWheels = Number.isFinite(params.nDrivenWheels)
    ? Math.max(1, Math.floor(params.nDrivenWheels))
    : 8;

  const mEff = mTrain / nDrivenWheels;
  const I_eff = I_disc + mEff * R * R;

  const FdriveMax = muAdh * mEff * g;
  const tauDriveMax = FdriveMax * R;

  // Startbedingungen
  const hasV0 = Number.isFinite(params.v0);
  const v0 = hasV0 ? Math.max(0, params.v0) : null;
  const wheelFrequency = Number.isFinite(params.wheelFrequency)
    ? Math.max(0, params.wheelFrequency)
    : 12;

  let omega = hasV0 ? v0 / R : 2 * Math.PI * wheelFrequency;
  if (scenario === "Acceleration with Brake" && !hasV0) omega = 0;
  let v_mps = omega * R;

  // Bremsen (Pad-Rad)
  const nBK = Number.isFinite(params.nBK) ? Math.max(0, Math.floor(params.nBK)) : 1;
  const bkAngles = Array.isArray(params.bkAngles) ? params.bkAngles : [0];

  const mu = Number.isFinite(params.mu) ? Math.max(0, params.mu) : 0.4;
  const FN = Number.isFinite(params.FN) ? Math.max(0, params.FN) : 3000;
  const mu_s = Number.isFinite(params.mu_s) ? Math.max(mu, params.mu_s) : mu;

  const frictionToHeat = Number.isFinite(params.frictionToHeat)
    ? Math.max(0, Math.min(1, params.frictionToHeat))
    : 0.9;

  // Zusatzenergie pro Fläche (J/m²)
  const Edensity = Number.isFinite(params.Edensity) ? Math.max(0, params.Edensity) : 0;

  // Singularitätsfix: Quellradius ~ L/4
  const sourceRadius = Number.isFinite(params.sourceRadius)
    ? Math.max(0, params.sourceRadius)
    : Math.max(1e-6, padLengthClamped / 4);

  const rEffSq = distSq + sourceRadius * sourceRadius;

  // 2-Körper-Split
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

  // Verlust (kLoss) als Kernel-Dämpfung
  const h_loss = Number.isFinite(params.kLoss) ? Math.max(0, params.kLoss) : 200;
  const Cw_local = Math.max(1e-12, rho_w * cp_w * V_effect);
  const lambdaLoss = (h_loss * contactArea) / Cw_local;

  // Halbraum-Faktor (Oberfläche)
  const surfaceFactor = 2;

  // Buffers
  const pulseWheelSeq = new Float64Array(steps);
  const pulseRailSeq  = new Float64Array(steps);
  const velocityProfile = new Float64Array(steps);

  let theta = 0;
  const omegaEps = 1e-6;

  const NsubBase = 8;
  const tauFromAccel = (aCmd) => I_eff * (aCmd / R);

  for (let i = 0; i < steps; i++) {
    const dThetaStepEst = Math.abs(omega * dt);
    const Nsub = Math.max(
      NsubBase,
      Math.ceil(dThetaStepEst / Math.max(1e-12, contactAngle * 0.25))
    );
    const dtSub = dt / Nsub;

    for (let s = 0; s < Nsub; s++) {
      let tauBrakeKinetic = 0;
      let tauBrakeStaticMax = 0;
      let inAnyContact = false;

      const thetaWrapped = wrapAngle(theta);

      for (let j = 0; j < nBK; j++) {
        const bAngle = Number.isFinite(bkAngles[j]) ? bkAngles[j] : 0;
        const d = wrapAngle(thetaWrapped - wrapAngle(bAngle));
        const inContact = Math.abs(d) <= contactAngle * 0.5;
        if (!inContact) continue;

        inAnyContact = true;
        tauBrakeKinetic += mu * FN * R;
        tauBrakeStaticMax += mu_s * FN * R;
      }

      // τ_drive vor Adhäsionslimit
      let tauDriveTarget = 0;

      if (scenario === "Constant Speed") {
        tauDriveTarget = inAnyContact ? tauBrakeKinetic : 0;
      } else if (scenario === "Emergency Brake") {
        tauDriveTarget = tauFromAccel(3 * aNeg); // aNeg < 0
      } else if (scenario === "Braking") {
        tauDriveTarget = tauFromAccel(aNeg);
      } else if (scenario === "Acceleration") {
        tauDriveTarget = tauFromAccel(aPos);
      } else if (scenario === "Acceleration with Brake") {
        tauDriveTarget = tauFromAccel(aPos) + (inAnyContact ? tauBrakeKinetic : 0);
      }

      // Adhäsionslimit: nur positives Antriebsmoment begrenzen
      let tauDrive = tauDriveTarget;
      if (tauDrive > tauDriveMax) tauDrive = tauDriveMax;

      // Haftreibung (nicht für "Acceleration with Brake")
      if (
        scenario !== "Acceleration with Brake" &&
        inAnyContact &&
        Math.abs(omega) < omegaEps
      ) {
        const tauDrivePos = Math.max(0, tauDrive);
        if (tauDrivePos <= tauBrakeStaticMax) {
          omega = 0;
          v_mps = 0;
          continue;
        }
      }

      const tauBrake = inAnyContact ? tauBrakeKinetic : 0;
      const tauNet = tauDrive - tauBrake;
      const alphaRot = tauNet / I_eff;

      theta = theta + omega * dtSub + 0.5 * alphaRot * dtSub * dtSub;

      omega = omega + alphaRot * dtSub;
      if (omega < 0) omega = 0;
      v_mps = omega * R;

      if (inAnyContact && omega > omegaEps) {
        const P_fric = tauBrake * omega;
        const dE_heat_total = frictionToHeat * P_fric * dtSub;

        const slipDist = v_mps * dtSub;
        const frac = Math.max(0, Math.min(1, slipDist / contactArc));
        const dE_absorb_total = Edensity * contactArea * frac;

        const dE_total = dE_heat_total + dE_absorb_total;

        pulseWheelSeq[i] += heatSplitWheel * dE_total;
        pulseRailSeq[i]  += heatSplitRail  * dE_total;
      }
    }

    velocityProfile[i] = v_mps;
  }

  const powerSeq = new Float64Array(steps);
  for (let i = 0; i < steps; i++) powerSeq[i] = pulseWheelSeq[i] / dt;

  // Green-Kernel
  const kernel = new Float64Array(steps);
  for (let i = 0; i < steps; i++) {
    const tau = (i + 0.5) * dt;
    const denom = Math.pow(4 * Math.PI * alpha_w * tau, 1.5);
    const G =
      denom > 0
        ? Math.exp(-rEffSq / (4 * alpha_w * tau)) / (denom * rho_w * cp_w)
        : 0;

    kernel[i] = surfaceFactor * G * Math.exp(-lambdaLoss * tau);
  }

  const conv = convolveFFT_real(pulseWheelSeq, kernel);

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
      // (Debug-Infos kannst du hier bei Bedarf wieder ergänzen)
    });
  }

  return out;
}
