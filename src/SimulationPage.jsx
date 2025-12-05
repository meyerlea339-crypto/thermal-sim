import { useState, useEffect } from "react";
import { useSimulationWorker } from "./lib/useSimulationWorker";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer
} from "recharts";
import { Settings2, Thermometer, Zap, Gauge, BatteryCharging } from "lucide-react";

const DEFAULT_PARAMS = {
  scenario: "Emergency Brake",   // Bremsung testen
  dt: 0.001,                     // 1 ms Zeitschritt (gute Auflösung, noch schnell)
  duration: 20,                  // 20 Sekunden Gesamtzeit
  wheelFrequency: 12,            // Startdrehzahl ~ 12 Hz (~720 U/min)
  aNeg: -1,                      // kräftige Verzögerung [m/s²]
  aPos: 2,                       // (nur für Acceleration relevant)
  kLoss: 0.05,                    // moderater Wärmeverlustfaktor 1/s
  Edensity: 0.05,                 // nicht genutzt bei Bremsung, nur Constant Speed
  thickness: 0.03,                // Scheibendicke 3 cm
  enableFiniteThickness: true,    // Dickenkorrektur aktiv
  rho_disc: 7800,                 // Stahl [kg/m³]
  cp_disc: 460,                   // Wärmekapazität J/(kg*K)
  k_disc: 50,                     // Wärmeleitfähigkeit W/(mK)
  distance: 0.05,                 // Messpunkt 5 cm vom Bremsklotz
  nBK: 2,                         // zwei Bremsklötze
  bkAngles: [0, Math.PI],         // gegenüberliegende Positionen
  tFirst: 0.005,                  // Startzeit erstes Ereignis
  R: 0.3,                         // Rad-Radius [m]
  mu: 0.4,                        // Reibungskoeffizient
  FN: 500,                        // Normalkraft pro Bremsklotz (N)
  massWheel: 80,                  // Radmasse in kg
  maxPoints: 2000,                // Chart-Auflösung / Performance
  T0: 20                          // Umgebungstemperatur °C
};

export default function SimulationPage() {
  const [params, setParams] = useState(DEFAULT_PARAMS);
  const [showFlash, setShowFlash] = useState(false);

  const data = useSimulationWorker(params);

  // Automatische kLoss-Anpassung (optional)
  useEffect(() => {
    let newKLoss = params.kLoss;
    if (params.scenario === "Constant Speed") newKLoss = 0.1;
    if (params.scenario === "Emergency Brake") newKLoss = 0.5;
    if (params.scenario === "Acceleration with Brake") newKLoss = 0.3;
    setParams((p) => ({ ...p, kLoss: newKLoss }));
  }, [params.scenario]);

  const updateParam = (key, value) => {
    setParams(prev => ({ ...prev, [key]: value }));
  };

  return (
    <div className="flex h-screen w-full bg-slate-950 text-slate-200 overflow-hidden font-sans">
      {/* Sidebar für Parameter */}
      <aside className="w-64 flex-shrink-0 border-r border-slate-800 bg-slate-900 overflow-y-auto p-4">
        <div className="flex items-center gap-2 mb-6 text-blue-400">
          <Settings2 className="w-5 h-5" />
          <h1 className="font-bold text-lg">Thermal Lab</h1>
        </div>

        {/* Szenario */}
        <Section title="Scenario">
          <select
            className="w-full bg-slate-800 border border-slate-700 rounded px-3 py-2 text-sm"
            value={params.scenario}
            onChange={(e) => updateParam("scenario", e.target.value)}
          >
            <option value="Constant Speed">Constant Speed</option>
            <option value="Emergency Brake">Emergency Brake</option>
            <option value="Acceleration with Brake">Acceleration with Brake</option>
          </select>
        </Section>

        {/* Systemparameter */}
        <Section title="System Constants">
          <InputGroup label="Duration [s]" value={params.duration} onChange={(v) => updateParam("duration", v)} min={1} max={100} step={1} />
          <InputGroup label="dt [s]" value={params.dt} onChange={(v) => updateParam("dt", v)} min={0.0005} max={0.02} step={0.0005} />
          <InputGroup label="Dist [m]" value={params.distance} onChange={(v) => updateParam("distance", v)} min={0.0001} max={0.01} step={0.0001} />
          <InputGroup label="Thermal Diff [m²/s]" value={params.thermalDiffusion} onChange={(v) => updateParam("thermalDiffusion", v)} min={1e-6} max={1e-4} step={1e-7} scientific />
          <InputGroup label="Convection Factor" value={params.kLoss} onChange={(v) => updateParam("kLoss", v)} min={0} max={1} step={0.01} />
          <InputGroup label="Pulse Energy [J]" value={params.Edensity} onChange={(v) => updateParam("Edensity", v)} min={0} max={1} step={0.01} />
          <InputGroup label="Wheel Freq [Hz]" value={params.wheelFrequency} onChange={(v) => updateParam("wheelFrequency", v)} min={0} max={20} step={0.1} />
          <InputGroup label="Wheel Mass [kg]" value={params.massWheel} onChange={(v) => updateParam("massWheel", v)} min={1} max={500} step={1} />
        </Section>
      </aside>

      {/* Hauptbereich mit Charts */}
      <main className="flex-1 p-6 space-y-6 overflow-auto">
        <Chart title="Wheel Point Temperature [°C]" icon={<Thermometer />} dataKey="bulkTemp" stroke="#f97316"
          data={data} showFlash={showFlash} setShowFlash={setShowFlash} height={300} />
        <Chart title="Cumulative Input Energy [J]" icon={<BatteryCharging />} dataKey="cumulativeEnergy" stroke="#facc15"
          data={data} height={250} />
        <Chart title="Input Power [W]" icon={<Zap />} dataKey="power" stroke="#ef4444"
          data={data} type="step" height={250} />
        <Chart title="Velocity [m/s]" icon={<Gauge />} dataKey="velocity" stroke="#10b981"
          data={data} height={250} />
      </main>
    </div>
  );
}

function Section({ title, children }) {
  return (
    <div className="space-y-4 border-t border-slate-800 pt-4">
      <h3 className="text-xs font-semibold text-slate-500 uppercase">{title}</h3>
      {children}
    </div>
  );
}

function InputGroup({ label, value, onChange, min, max, step, scientific = false }) {
  const safeVal = Number.isFinite(value) ? value : 0;
  return (
    <div className="space-y-1">
      <div className="flex justify-between">
        <label className="text-xs text-slate-400">{label}</label>
        <span className="text-xs text-slate-500 font-mono">
          {scientific ? safeVal.toExponential(2) : safeVal.toFixed(3)}
        </span>
      </div>
      <div className="flex items-center gap-2">
        <input type="range" min={min} max={max} step={step}
          value={safeVal} onChange={(e) => onChange(parseFloat(e.target.value))}
          className="flex-1 h-1.5 bg-slate-700 rounded-lg accent-blue-500" />
        <input type="number" value={safeVal} step={step}
          onChange={(e) => onChange(parseFloat(e.target.value))}
          className="w-16 bg-slate-800 border border-slate-700 rounded px-1 py-0.5 text-xs text-right focus:border-blue-500" />
      </div>
    </div>
  );
}

function Chart({ title, icon, dataKey, stroke, data, showFlash, setShowFlash, type = "monotone", height }) {
  const safeData = Array.isArray(data) ? data : [];
  const xDomain = safeData.length ? [0, Math.max(...safeData.map(d => d.time || 0))] : [0, 1];
  const yDomain = safeData.length ? ['dataMin', 'dataMax'] : [0, 1];

  return (
    <div className="bg-slate-900 border border-slate-800 rounded-lg p-4 shadow-sm w-full mb-6">
      <div className="flex items-center justify-between mb-4 text-orange-400">
        <div className="flex items-center gap-2">{icon} <h3 className="font-medium text-sm">{title}</h3></div>
        {setShowFlash && (
          <div className="flex items-center gap-2 text-xs text-slate-400">
            <label>Show Flash Pulses</label>
            <input type="checkbox" checked={showFlash}
              onChange={(e) => setShowFlash(e.target.checked)} />
          </div>
        )}
      </div>
      <ResponsiveContainer width="100%" height={height} minWidth={0}>
        <LineChart data={safeData}>
          <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
          <XAxis type="number" dataKey="time" domain={xDomain} stroke="#94a3b8" fontSize={12} />
          <YAxis domain={yDomain} stroke="#94a3b8" fontSize={12} />
          <Tooltip />
          <Legend />
          <Line type={type} dataKey={dataKey} stroke={stroke} strokeWidth={2} dot={false} />
          {showFlash && <Line type="monotone" dataKey="combinedTemp"
            stroke="#ef4444" strokeWidth={1} strokeOpacity={0.5} dot={false} />}
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
}