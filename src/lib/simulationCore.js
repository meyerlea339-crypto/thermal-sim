// ----- FFT-basierte schnelle Faltung ---------------------------------------

const safe = (v, fb = 0) => Number.isFinite(v) ? v : fb;

// In-place Radix-2 FFT
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
      for (let j2 = 0; j2 < len / 2; j2++) {
        const uRe = re[i + j2];
        const uIm = im[i + j2];

        const vRe = re[i + j2 + len/2] * wRe - im[i + j2 + len/2] * wIm;
        const vIm = re[i + j2 + len/2] * wIm + im[i + j2 + len/2] * wRe;

        re[i + j2]             = uRe + vRe;
        im[i + j2]             = uIm + vIm;
        re[i + j2 + len/2] = uRe - vRe;
        im[i + j2 + len/2] = uIm - vIm;

        const nwRe = wRe * wlenRe - wIm * wlenIm;
        const nwIm = wRe * wlenIm + wIm * wlenRe;
        wRe = nwRe; wIm = nwIm;
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
  const na = a.length, nb = b.length;
  const outLen = na + nb - 1;

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

// --------------------------------------------------------------------------
// ---------   F F T   S I M U L A T I O N   C O R E   ----------------------
// --------------------------------------------------------------------------

export function simulationCore(params) {
  const dt = Math.max(1e-4, params.dt || 1e-3);
  const duration = Math.max(dt, params.duration || 5);
  const steps = Math.floor(duration / dt);

  // --- Material & Geometrie ---
  const rho = params.rho_disc;
  const cp = params.cp_disc;
  const k = params.k_disc;
  const R = params.R;
  const thickness = params.thickness;
  const volume = Math.PI * R * R * thickness;
  const mass = volume * rho;
  const alpha = k / (rho * cp);
  const distSq = params.distance * params.distance;
  const T0 = params.T0 || 20;

  // --- Wärmeprofil für 1 Joule (Impulsantwort) ---
  const pulseProfile = new Float64Array(steps);
  for (let i = 0; i < steps; i++) {
    const tau = i * dt;
    if (tau > 0) {
      const temporalTerm = Math.pow(4 * Math.PI * alpha * tau, 1.5);
      const spatialTerm = Math.exp(-distSq / (4 * alpha * tau));
      pulseProfile[i] = (1 / (mass * cp)) * spatialTerm / temporalTerm;
    }
  }

  // --- Initialisierung ---
  const pulseSeq = new Float64Array(steps);
  const velocityProfile = new Float64Array(steps); // m/s

  let freq = params.scenario === "Acceleration with Brake"
    ? 0 // Start bei Stillstand
    : (params.wheelFrequency || 10); // Hertz für andere Szenarien

  const FN = params.FN;
  const mu = params.mu;
  const aNeg = params.aNeg || -2;   // m/s² (Bremsverzögerung)
  const aPos = params.aPos || 1;    // m/s² (Antrieb)
  const omegaConst = 2 * Math.PI;
  const nBK = params.nBK;
  const anglePerRev = 2 * Math.PI;
  let cumAng = 0;
  const nextPulse = params.bkAngles.slice();

  // --- Haupt-Simulationsschleife ---
  for (let i = 0; i < steps; i++) {
    const v_mps = freq * (2 * Math.PI * R); // aktuelle lineare Geschwindigkeit [m/s]

    if (params.scenario === "Emergency Brake") {
      // reine Bremsung
      const new_v = Math.max(v_mps + aNeg * dt, 0);
      freq = new_v / (2 * Math.PI * R);

    } else if (params.scenario === "Acceleration with Brake") {
      // Beschleunigung + Bremse gleichzeitig
      const net_a = aPos + aNeg; // aNeg ist negativ!
      const new_v = Math.max(v_mps + net_a * dt, 0);
      freq = new_v / (2 * Math.PI * R);

    } else if (params.scenario === "Constant Speed") {
      // Geschwindigkeit bleibt konstant
      // freq unverändert
    }

    // Speichern der linearen Geschwindigkeit in m/s
    velocityProfile[i] = freq * (2 * Math.PI * R);

    // --- Pulsereignisse berechnen ---
    cumAng += freq * omegaConst * dt;
    for (let j = 0; j < nBK; j++) {
      if (cumAng >= nextPulse[j]) {
        const M = mu * FN * R;              // Reibungsmoment [Nm]
        const Q = M * (Math.PI / nBK);      // Energie pro Puls [J]
        pulseSeq[i] += Q;
        nextPulse[j] += anglePerRev;        // nächster Puls
      }
    }
  }

  // --- Temperaturantwort berechnen ---
  const tempResponse = new Float64Array(steps);
  for (let t = 0; t < steps; t++) {
    let sum = 0;
    for (let p = 0; p <= t; p++) {
      const Q = pulseSeq[p];
      if (Q > 0) {
        sum += Q * pulseProfile[t - p] * Math.exp(-params.kLoss * (t - p) * dt);
      }
    }
    tempResponse[t] = T0 + sum;
  }

  // --- Ergebnis zusammenstellen ---
  const samplingStride = Math.max(1, Math.floor(steps / (params.maxPoints || 2000)));
  const result = [];
  for (let i = 0; i < steps; i += samplingStride) {
    result.push({
      time: i * dt,
      bulkTemp: tempResponse[i],
      cumulativeEnergy: pulseSeq.slice(0, i + 1).reduce((a, b) => a + b, 0),
      power: pulseSeq[i] / dt,  // Instantanleistung [W]
      velocity: velocityProfile[i] // m/s
    });
  }

  return result;
}