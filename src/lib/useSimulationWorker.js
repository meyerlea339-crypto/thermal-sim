import { useEffect, useRef, useState } from "react";

// Custom-Hook für Web Worker
export function useSimulationWorker(params) {
  const [data, setData] = useState([]);
  const workerRef = useRef(null);

  useEffect(() => {
    // Worker starten
    workerRef.current = new Worker(
      new URL("./simulation.worker.js", import.meta.url),
      { type: "module" }
    );

    // Nachricht vom Worker empfangen
    workerRef.current.onmessage = (event) => {
      const msg = event.data;
      // Nur Arrays annehmen (fehlerfrei für Recharts)
      if (Array.isArray(msg)) {
        setData(msg);
      } else {
        setData([]);
      }
    };

    // Worker beim Unmount beenden
    return () => {
      workerRef.current.terminate();
    };
  }, []);

  // Simulation starten bei Param-Änderung
  useEffect(() => {
    if (workerRef.current) {
      workerRef.current.postMessage(params);
    }
  }, [params]);

  return Array.isArray(data) ? data : [];
}