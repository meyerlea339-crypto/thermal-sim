import { useEffect, useRef, useState } from "react";

const safeNum = (v, fb = 0) => {
  const n = Number(v);
  return Number.isFinite(n) ? n : fb;
};

export function useSimulationWorker(params) {
  const [data, setData] = useState([]);
  const workerRef = useRef(null);

  useEffect(() => {
    workerRef.current = new Worker(
      new URL("./simulation.worker.js", import.meta.url),
      { type: "module" }
    );

    workerRef.current.onmessage = (event) => {
      const msg = event.data;
      if (!Array.isArray(msg)) {
        setData([]);
        return;
      }
      const sanitized = msg.map((d) => ({
        time: safeNum(d.time, 0),
        pointTemp: safeNum(d.pointTemp, 0), // wichtig: PointTemp aus simulationCore
        cumulativeEnergy: safeNum(d.cumulativeEnergy, 0),
        power: safeNum(d.power, 0),
        velocity: safeNum(d.velocity, 0)
      }));
      setData(sanitized);
    };

    return () => {
      workerRef.current && workerRef.current.terminate();
      workerRef.current = null;
    };
  }, []);

  useEffect(() => {
    if (workerRef.current) {
      workerRef.current.postMessage(params);
    }
  }, [params]);

  return Array.isArray(data) ? data : [];
}