// Schnelle Greensche Wärmeleitung für 3D-Punktquelle
// Optimiert auf O(N) Laufzeit durch akkumulierten Temperaturenergie-Pool

export function useSimulation(params) {
  const {
    scenario,
    thermalDiffusion, // α [m²/s] - kann auch als Effektparameter genutzt werden
    distance, // r [m] - Abstand vom Wärmeeintragpunkt
    dt,
    duration,
    wheelFrequency,
    aNeg,
    aPos,
    Edensity,
    kLoss,
    thickness,
    enableFiniteThickness,
    rho_disc,
    cp_disc,
    k_disc
  } = params;

  const data = [];

  const T0 = 20; // Umgebungstemperatur in °C
  let velocity = wheelFrequency; // aktuelle Geschwindigkeit
  let cumulativeEnergy = 0;

  // Thermische Diffusivität aus Materialparametern
  const alpha_disc = k_disc / (rho_disc * cp_disc);

  // Akkumulierter Temperaturenergie-Pool - speichert Energieeintrag + Zerfall
  let tempEnergyPool = 0;

  for (let t = 0; t <= duration; t += dt) {
    let power = 0;

    // === Szenarien-Formeln ===
    if (scenario === "Constant Speed") {
      // konstante Energiezufuhr
      power = Edensity * wheelFrequency;
    } else if (scenario === "Emergency Brake") {
      velocity += aNeg * dt;
      if (velocity < 0) velocity = 0;
      power = Math.abs(aNeg) * 15; // mehr Energie durch starkes Bremsen
    } else if (scenario === "Acceleration with Brake") {
      velocity += aPos * dt;
      power = aPos * 5; // Energieeintrag durch Beschleunigen
      if (t % 2 < dt) {
        power += 40; // kleiner Bremsimpuls alle 2 Sekunden
      }
    }

    // === Energie dieses Schritts ===
    const Q_step = power * dt;
    cumulativeEnergy += Q_step;

    // === Temperaturzuwachs durch aktuellen Puls ===
    // Green’sche Funktion für 3D-Punktquelle:
    // ΔT = Q / (ρ c_p (4 π α τ)^(3/2)) * exp(- r² / (4 α τ))
    const deltaTemp =
      Q_step /
      (rho_disc * cp_disc * Math.pow(4 * Math.PI * alpha_disc * dt, 1.5)) *
      Math.exp(-(distance ** 2) / (4 * alpha_disc * dt));

    // === Akkumulieren und Zerfall anwenden ===
    tempEnergyPool = tempEnergyPool * Math.exp(-kLoss * dt) + deltaTemp;

    // === Gesamttemperatur ===
    let bulkTemp = T0 + tempEnergyPool;

    // Option: Finite Thickness beeinflusst Wärmeverlust durch Diffusion
    if (enableFiniteThickness) {
      const diffusionFactor = thermalDiffusion * dt / thickness;
      bulkTemp -= diffusionFactor * (bulkTemp - T0);
    }

    // Oberflächentemperatur inkl. "Flash" Oszillation
    const combinedTemp = bulkTemp + Math.sin(t * 8) * 1.5;

    // === Speichern der aktuellen Werte ===
    data.push({
      time: parseFloat(t.toFixed(3)),
      bulkTemp,
      combinedTemp,
      cumulativeEnergy,
      power,
      velocity
    });
  }

  return data;
}