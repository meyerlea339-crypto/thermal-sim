// SimulationPage.jsx
import { useState, useLayoutEffect } from "react";
import { useSimulationWorker } from "./lib/useSimulationWorker";
import {
  ResponsiveContainer,
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend
} from "recharts";
import {
  Settings2,
  Thermometer,
  Zap,
  Gauge,
  BatteryCharging
} from "lucide-react";

const DEFAULT_PARAMS = {
  scenario: "Constant Speed",
  dt: 0.001,
  duration: 15,
  wheelFrequency: 12,
  aNeg: -2.5,
  aPos: 2.0,
  Edensity: 50,
  thermalVolumeFactor: 1.0,
  k_disc: 50,
  rho_disc: 7800,
  cp_disc: 460,
  kLoss: 50,
  distance: 0.001,
  nBK: 1,
  bkAngles: [0, Math.PI],
  R: 0.3,
  mu: 0.4,
  FN: 3000,
  frictionToHeat: 0.9,
  wheelMass: 200,
  maxPoints: 500,
  T0: 20
};

/* Bildschirmgröße für spätere UI-Anpassungen (falls benötigt) */
function useDisplaySize() {
  const [size, setSize] = useState({ width: 0, height: 0 });

  useLayoutEffect(() => {
    if (typeof window === "undefined") return;

    const measure = () => {
      setSize({
        width: window.innerWidth,
        height: window.innerHeight
      });
    };

    measure();
    window.addEventListener("resize", measure);
    return () => window.removeEventListener("resize", measure);
  }, []);

  return size;
}

export default function SimulationPage() {
  const [params, setParams] = useState(DEFAULT_PARAMS);
  const data = useSimulationWorker(params);

  const updateParam = (key, value) =>
    setParams(prev => ({ ...prev, [key]: value }));

  // velocity in km/h
  const chartData = Array.isArray(data)
    ? data.map(d => ({
        ...d,
        velocity:
          Number.isFinite(d.velocity) ? d.velocity * 3.6 : d.velocity
      }))
    : [];

  return (
    <div className="flex h-screen w-screen bg-slate-950 text-slate-200 overflow-hidden font-sans">
      {/* Sidebar jetzt breiter */}
      <aside className="w-80 flex-none border-r border-slate-800 bg-slate-900 overflow-y-auto p-4">
        <div className="flex items-center gap-2 mb-6 text-blue-400">
          <Settings2 className="w-6 h-6" />
          <h1 className="font-bold text-2xl">Thermal Lab</h1>
        </div>

        <Section title="Scenario">
          <select
            value={params.scenario}
            onChange={(e) => updateParam("scenario", e.target.value)}
            className="w-full bg-slate-800 border border-slate-700 rounded px-3 py-2 text-sm"
          >
            <option value="Constant Speed">Constant Speed</option>
            <option value="Emergency Brake">Emergency Brake</option>
            <option value="Acceleration with Brake">Acceleration with Brake</option>
          </select>
        </Section>

        <Section title="System Constants">
          <InputGroup label="Duration [s]" value={params.duration} onChange={(v) => updateParam("duration", v)} min={1} max={120} step={0.5} />
          <InputGroup label="dt [s]" value={params.dt} onChange={(v) => updateParam("dt", v)} min={1e-6} max={0.01} step={1e-6} scientific />
          <InputGroup label="Distance [m]" value={params.distance} onChange={(v) => updateParam("distance", v)} min={0} max={0.02} step={0.0001} />
          <InputGroup label="Thermal k [W/mK]" value={params.k_disc} onChange={(v) => updateParam("k_disc", v)} min={1} max={200} step={1} />
          <InputGroup label="Energie pro Kontaktfläche [J/m²]" value={params.Edensity} onChange={(v) => updateParam("Edensity", v)} min={0} max={1000} step={1} />
          <InputGroup label="Wheel Freq [Hz]" value={params.wheelFrequency} onChange={(v) => updateParam("wheelFrequency", v)} min={0} max={50} step={0.5} />
          <InputGroup label="mu" value={params.mu} onChange={(v) => updateParam("mu", v)} min={0} max={2} step={0.01} />
          <InputGroup label="FN [N]" value={params.FN} onChange={(v) => updateParam("FN", v)} min={0} max={50000} step={10} />
          <InputGroup label="Bremsklötze" value={params.nBK} onChange={(v) => updateParam("nBK", Math.max(1, Math.min(2, v)))} min={1} max={2} step={1} />
          <InputGroup label="Wärmeübergang h [W/m²K]" value={params.kLoss} onChange={(v) => updateParam("kLoss", v)} min={0} max={2000} step={1} />
        </Section>
      </aside>

      {/* Chartbereich — füllt den Rest komplett aus */}
      <main className="flex-1 p-4 overflow-hidden min-h-0 min-w-0 h-full">
        <div className="grid grid-cols-2 grid-rows-2 gap-6 h-full">
          <Chart
            title="Wheel Point Temperature [°C]"
            icon={<Thermometer />}
            dataKey="pointTemp"
            stroke="#f97316"
            data={chartData}
            clampLowerToZero
          />
          <Chart
            title="Cumulative Input Energy [J]"
            icon={<BatteryCharging />}
            dataKey="cumulativeEnergy"
            stroke="#facc15"
            data={chartData}
            clampLowerToZero
          />
          <Chart
            title="Input Power [W]"
            icon={<Zap />}
            dataKey="power"
            stroke="#ef4444"
            data={chartData}
            type="step"
            clampLowerToZero
          />
          <Chart
            title="Velocity [km/h]"
            icon={<Gauge />}
            dataKey="velocity"
            stroke="#10b981"
            data={chartData}
            domainAuto
          />
        </div>
      </main>
    </div>
  );
}

/* ---------- UI Hilfs-Komponenten ---------- */

function Section({ title, children }) {
  return (
    <div className="space-y-4 border-t border-slate-800 pt-4 mb-6">
      <h3 className="text-xs font-semibold text-slate-500 uppercase tracking-wide">
        {title}
      </h3>
      {children}
    </div>
  );
}

function InputGroup({
  label,
  value,
  onChange,
  min,
  max,
  step,
  scientific = false
}) {
  const safe = Number.isFinite(value) ? value : 0;

  return (
    <div className="space-y-1">
      <div className="flex justify-between items-center">
        <label className="text-xs text-slate-400">{label}</label>
        <span className="text-xs text-slate-500 font-mono">
          {scientific ? safe.toExponential(2) : safe.toFixed(3)}
        </span>
      </div>

      <div className="flex items-center gap-2">
        <input
          type="range"
          min={min}
          max={max}
          step={step}
          value={safe}
          onChange={(e) => onChange(parseFloat(e.target.value))}
          className="flex-1 h-1.5 bg-slate-700 rounded-lg accent-blue-500"
        />
        <input
          type="number"
          value={safe}
          step={step}
          onChange={(e) => onChange(parseFloat(e.target.value))}
          className="w-20 bg-slate-800 border border-slate-700 rounded px-1 py-0.5 text-xs text-right"
        />
      </div>
    </div>
  );
}

/* ---------- Chart-Komponente ---------- */

function Chart({
  title,
  icon,
  dataKey,
  stroke,
  data,
  type = "monotone",
  clampLowerToZero = false,
  domainAuto = false
}) {
  // Sortieren nach time
  const safeData = Array.isArray(data)
    ? [...data].sort((a, b) => (a.time || 0) - (b.time || 0))
    : [];

  const xMax = safeData.length
    ? Math.max(...safeData.map((d) => d.time || 0))
    : 1;
  const xDomain = [0, xMax];

  // Y-Achse
  const vals = safeData
    .map((d) => Number(d[dataKey]))
    .filter(Number.isFinite);

  let yMin = vals.length ? Math.min(...vals) : 0;
  let yMax = vals.length ? Math.max(...vals) : 1;

  if (yMin === yMax) {
    const eps = Math.abs(yMin) * 0.02 || 1;
    yMin -= eps;
    yMax += eps;
  }

  const pad = (yMax - yMin) * 0.1;
  let lower = yMin - pad;
  let upper = yMax + pad;

  if (clampLowerToZero) lower = Math.max(0, lower);
  if (upper <= lower) upper = lower + Math.abs(lower) * 0.02;

  const yDomain = [lower, upper];

  return (
    <div className="bg-slate-900 border border-slate-800 rounded-lg overflow-hidden flex flex-col w-full h-full">
      <div className="flex items-center justify-between px-4 py-3 text-orange-400">
        <div className="flex items-center gap-2">
          {icon}
          <h3 className="font-medium text-sm">{title}</h3>
        </div>
      </div>

      <div className="flex-1 px-2 pb-3 min-h-0 min-w-0">
        <ResponsiveContainer width="100%" height="100%">
          <LineChart
            data={safeData}
            margin={{ top: 8, right: 12, left: 0, bottom: 36 }}
          >
            <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
            <XAxis
              type="number"
              dataKey="time"
              domain={xDomain}
              stroke="#94a3b8"
              fontSize={12}
              tickFormatter={(v) => v.toFixed(2)}
              height={28}
              tickMargin={8}
            />
            <YAxis
              domain={domainAuto ? ["auto", "auto"] : yDomain}
              stroke="#94a3b8"
              fontSize={12}
              tickFormatter={(v) =>
                Math.abs(v) >= 1000 ? (v / 1000).toFixed(1) + "k" : v.toFixed(2)
              }
            />
            <Tooltip
              formatter={(v) =>
                Number.isFinite(v) ? Number(v.toFixed(6)) : v
              }
            />
            <Legend
              verticalAlign="top"
              align="right"
              height={28}
            />
            <Line
              type={type}
              dataKey={dataKey}
              stroke={stroke}
              strokeWidth={2}
              dot={false}
              connectNulls
              isAnimationActive={false}
            />
          </LineChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
