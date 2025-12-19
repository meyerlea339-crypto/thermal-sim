// simulation.worker.js
import { simulationCore } from "./simulationCore";

self.onmessage = (e) => {
  try {
    // Hook schickt direkt die params:
    const params = e.data || {};
    const result = simulationCore(params);
    // Hook erwartet direkt ein Array als Nachricht:
    self.postMessage(result);
  } catch (err) {
    console.error("Simulation worker error:", err);
    // Im Fehlerfall leeres Array senden:
    self.postMessage([]);
  }
};
