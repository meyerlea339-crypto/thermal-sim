/* eslint-disable no-restricted-globals */
import { simulationCore } from "./simulationCore.js";

// Worker hört auf Nachrichten vom Hauptthread
self.onmessage = (event) => {
  const params = event.data;
  try {
    // Simulation ausführen
    const result = simulationCore(params);
    // Ergebnis zurück an Hauptthread
    postMessage(Array.isArray(result) ? result : []);
  } catch (err) {
    console.error("Worker simulation error:", err);
    // Bei Fehler leeres Array senden
    postMessage([]);
  }
};