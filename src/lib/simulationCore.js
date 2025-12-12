// simulationCore.js
// Rotational wheel model — Edensity als J/m^2 (pro Kontaktfläche), realistischere Kühlung mittels h [W/m^2/K].
// Export: simulationCore(params) => Array of samples

function fft_inplace(re, im, invert = false) {
  const n = re.length;
  let j = 0;
  for (let i = 1; i < n; i++) {
    let bit = n >> 1;
    while (j & bit) { j ^= bit; bit >>= 1; }
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = 2 * Math.PI / len * (invert ? -1 : 1);
    const wlenRe = Math.cos(ang);
    const wlenIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wRe = 1, wIm = 0;
      const half = len >> 1;
      for (let j2 = 0; j2 < half; j2++) {
        const uRe = re[i + j2];
        const uIm = im[i + j2];
        const vr = re[i + j2 + half];
        const vi = im[i + j2 + half];
        const vRe = vr * wRe - vi * wIm;
        const vIm = vr * wIm + vi * wRe;
        re[i + j2] = uRe + vRe;
        im[i + j2] = uIm + vIm;
        re[i + j2 + half] = uRe - vRe;
        im[i + j2 + half] = uIm - vIm;
        const nwRe = wRe * wlenRe - wIm * wlenIm;
        const nwIm = wRe * wlenIm + wIm * wlenRe;
        wRe = nwRe; wIm = nwIm;
      }
    }
  }
  if (invert) {
    for (let i = 0; i < n; i++) {
      re[i] /= n; im[i] /= n;
    }
  }
}

function convolveFFT_real(a, b) {
  const na = a.length, nb = b.length;
  const outLen = na + nb - 1;
  let n = 1; while (n < outLen) n <<= 1;
  const are = new Float64Array(n);
  const aim = new Float64Array(n);
  const bre = new Float64Array(n);
  const bim = new Float64Array(n);
  are.set(a); bre.set(b);
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

export function simulationCore(params) {
  const dt = Math.max(1e-6, params.dt || 1e-3);
  const duration = Math.max(dt, params.duration || 5);
  const steps = Math.max(1, Math.floor(duration / dt));

  // Materialkonstanten
  const rho = Number.isFinite(params.rho_disc) ? params.rho_disc : 7800; // kg/m3
  const cp = Number.isFinite(params.cp_disc) ? params.cp_disc : 460;    // J/(kg K)
  let k = Number.isFinite(params.k_disc) ? params.k_disc : 50;          // W/(m K)
  if (k <= 0) k = 1e-12;
  const alpha = k / (rho * cp);
  const T0 = Number.isFinite(params.T0) ? params.T0 : 20; // °C

  // Rad / Szenario
  const R = Number.isFinite(params.R) ? params.R : 0.3; // m
  const scenario = params.scenario || "Constant Speed";
  const initialFreq = Number.isFinite(params.wheelFrequency) ? params.wheelFrequency : 10; // Hz
  const aNeg = Number.isFinite(params.aNeg) ? params.aNeg : -1; // m/s^2
  const aPos = Number.isFinite(params.aPos) ? params.aPos : 2;  // m/s^2
  const v_ref = Math.max(initialFreq, 1e-6) * 2 * Math.PI * R;

  // Kontakt-Geometrie
  const contactAngle = 0.05;         // rad, Bremsklotzbreite in Rad
  const contactDepth = 0.001;        // m, Eindringtiefe / Wärmetiefe (typisch sehr klein)
  const contactArea = contactAngle * R * 1.0; // m^2 (1 m axial angenommen)

  // Edensity: jetzt INTERPRETIERT ALS J/m^2 (Energie pro Kontaktfläche)
  // Benutzer-Parameter: params.Edensity [J/m^2]
  const Edensity = Number.isFinite(params.Edensity) ? params.Edensity : 0.05; // J/m^2 default (sehr klein)
  const E_absorb_base = Edensity * contactArea; // J per contact (bei speedScale = 1)

  // Effektives Volumen, das die Wärme initial aufnimmt (für Kühlung/Abkühlung)
  // Standard: contactVolume = contactArea * contactDepth
  const contactVolume = contactArea * contactDepth; // m^3
  // Option: allow scaling the effective thermal volume (user can increase to model deeper spread)
  const thermalVolumeFactor = Number.isFinite(params.thermalVolumeFactor) ? Math.max(1e-6, params.thermalVolumeFactor) : 1.0;
  const V_effect = contactVolume * thermalVolumeFactor; // m^3

  // Reib- / Brems-Parameter (mechanische Arbeit)
  const mu_param = Number.isFinite(params.mu) ? Math.max(0, params.mu) : 0;
  const FN_param = Number.isFinite(params.FN) ? Math.max(0, params.FN) : 0;
  const slidingDistance = R * contactAngle; // m
  const F_friction = mu_param * FN_param; // N
  const E_mech_per_contact = F_friction * slidingDistance; // J per contact
  const frictionToHeat = Number.isFinite(params.frictionToHeat) ? Math.max(0, Math.min(1, params.frictionToHeat)) : 1.0;

  // Kühlkoeffizient: kLoss interpretieren als Wärmeübergangskoeffizient h [W/(m^2 K)]
  // Realistische Werte: natürliche Konvektion ~ 5..25, erzwungene Konvektion 10..200, Kontakt/Kühlkörper deutlich größer.
  const h = Number.isFinite(params.kLoss) ? Math.max(0, params.kLoss) : 20; // W / (m^2 K)

  // Rad-Trägheitsmoment: I_disc oder compute from wheelMass (solid disc)
  let I_disc = Number.isFinite(params.I_disc) ? params.I_disc : null;
  if (!I_disc) {
    const wheelMass = Number.isFinite(params.wheelMass) ? params.wheelMass : 100; // kg
    I_disc = 0.5 * wheelMass * R * R;
  }

  // Pulse timing
  const nBK = Math.max(1, Math.min(2, Number.isFinite(params.nBK) ? params.nBK : 1));
  let bkAngles = nBK === 1 ? [0] : [0, Math.PI];
  const phi_obs = Number.isFinite(params.phi_obs) ? params.phi_obs : 0;
  bkAngles = bkAngles.map(a => ((a - phi_obs) + 2 * Math.PI) % (2 * Math.PI));
  const nextPulse = bkAngles.map(a => a);

  const dist = Number.isFinite(params.distance) ? params.distance : 0.001;
  const distSq = dist * dist;

  // Buffers
  const pulseSeq = new Float64Array(steps);     // Leistung W pro timestep, wird an Kernel gefaltet
  const velocityProfile = new Float64Array(steps);

  // Startbedingungen
  let freq = (scenario === "Acceleration with Brake") ? 0 : initialFreq;
  let omega = 2 * Math.PI * freq;
  let cumAngle = 0;

  // Simulation loop
  for (let i = 0; i < steps; i++) {
    // Szenario: kontinuierliche Beschleunigung/Verzögerung
    if (scenario === "Emergency Brake") {
      const alpha_cont = aNeg / R; // rad/s^2 (neg)
      omega = Math.max(0, omega + alpha_cont * dt);
      freq = omega / (2 * Math.PI);
    } else if (scenario === "Acceleration with Brake") {
      const alpha_drive = aPos / R;
      omega = Math.max(0, omega + alpha_drive * dt);
      freq = omega / (2 * Math.PI);
    } else {
      // Constant Speed: erzwinge exakt initialFreq (kein numerischer Drift)
      freq = initialFreq;
      omega = 2 * Math.PI * initialFreq;
    }

    const v_mps = omega * R;
    velocityProfile[i] = v_mps;

    const dAngle = freq * 2 * Math.PI * dt;
    const prevAngle = cumAngle;
    cumAngle += dAngle;

    // Kontakt-Events
    for (let j = 0; j < nBK; j++) {
      while (nextPulse[j] <= prevAngle) nextPulse[j] += 2 * Math.PI;
      if (prevAngle < nextPulse[j] && nextPulse[j] <= cumAngle) {
        // mechanische Bremsarbeit pro Kontakt (J)
        const E_brake_mech = E_mech_per_contact;

        // thermischer Anteil aus Reibung + Materialabsorption (Edensity per contact area)
        let speedScale = v_mps / (v_ref || 1e-9);
        if (!isFinite(speedScale) || speedScale < 0) speedScale = 0;
        const E_absorb = E_absorb_base * speedScale; // J per contact from material absorption (speed dependent)
        const E_heat_from_friction = frictionToHeat * E_brake_mech; // J from friction -> heat
        const E_pulse_thermal = E_absorb + E_heat_from_friction; // total J per contact into local thermal zone

        if (E_pulse_thermal > 0) {
          // convert to power for this timestep: W = J / dt
          pulseSeq[i] += E_pulse_thermal / dt;
        }

        // nur Energieentzug wenn nicht ConstantSpeed (Motor kompensiert sonst)
        if (scenario !== "Constant Speed") {
          const Ek_rot = 0.5 * I_disc * omega * omega;
          const Ek_new = Math.max(0, Ek_rot - E_brake_mech);
          omega = Ek_new > 0 ? Math.sqrt(2 * Ek_new / I_disc) : 0;
          freq = omega / (2 * Math.PI);
          velocityProfile[i] = omega * R;
        }

        nextPulse[j] += 2 * Math.PI;
      }
    }
  }

  // Green's kernel (simple 3D point-like kernel)
  const kernel = new Float64Array(steps);
  for (let i = 0; i < steps; i++) {
    const tau = (i === 0) ? dt * 0.5 : i * dt;
    const denom = Math.pow(4 * Math.PI * alpha * tau, 1.5);
    kernel[i] = (denom > 0)
      ? Math.exp(-distSq / (4 * alpha * tau)) / (denom * rho * cp)
      : 0;
  }

  // Convolution: pulseSeq (W) * kernel -> delta-T (K) time series
  let conv = convolveFFT_real(pulseSeq, kernel);

  // Kühlung: realistische Newton'sche Abkühlung auf das effektive Volumen V_effect
  // Wärmestrom q = h * A * (T - T0) -> Temperaturabfallrate dT/dt = - q / (m*cp) = - h*A/(rho*V*cp) * (T - T0)
  // => Lösung über dt: T(t+dt) = T0 + (T(t)-T0) * exp(-h*A/(rho*V*cp) * dt)
  const area_cool = contactArea; // Fläche für Wärmeabgabe (angepasst falls nötig)
  const coolingFactorPerSec = (h * area_cool) / (rho * V_effect * cp); // 1/s
  for (let i = 0; i < steps; i++) {
    let temp = T0 + conv[i]; // absolute T
    // exponential decay for stability
    temp = T0 + (temp - T0) * Math.exp(-coolingFactorPerSec * dt);
    conv[i] = temp - T0; // store delta-T
  }

  // ----- Stabiles Resampling für Ausgabe -----
  const maxPoints = Math.max(2, Math.min(5000, Number.isFinite(params.maxPoints) ? params.maxPoints : 500));
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
    const tFrac = (Nout === 1) ? 0 : (i / (Nout - 1));
    const tIndexFloat = tFrac * (steps - 1);
    const time = tFrac * duration;
    const pointTemp = T0 + sampleAt(conv, tIndexFloat);
    const cumulativeEnergy = (() => {
      const idx = Math.floor(tIndexFloat);
      let sum = 0;
      for (let j = 0; j <= idx; j++) sum += pulseSeq[j] * dt;
      const frac = tIndexFloat - idx;
      if (idx + 1 < pulseSeq.length) sum += pulseSeq[idx + 1] * dt * frac;
      return sum;
    })();
    const power = sampleAt(pulseSeq, tIndexFloat);
    const velocity = sampleAt(velocityProfile, tIndexFloat);

    out.push({
      time,
      pointTemp,
      cumulativeEnergy,
      power,
      velocity
    });
  }

  return out;
}