/* eslint-disable no-restricted-globals */
import { simulationCore } from "./simulationCore.js";

self.onmessage = (event) => {
  const params = event.data;

  try {
    const result = simulationCore(params);

    // Nur gültige Arrays an den Main Thread schicken
    postMessage(Array.isArray(result) ? result : []);
  } catch (err) {
    // Keine console.log hier nötig – Browser worker Konsole bleibt klein
    postMessage([]);
  }
};