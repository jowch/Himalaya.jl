import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { IndexEntry } from "../api";
import { phaseColor } from "../phases";

export interface MillerPlotProps {
  indices: IndexEntry[];
}

export interface ScatterRow { ratio: number; q: number; phase: string; color: string; }

export function toScatterData(indices: IndexEntry[]): ScatterRow[] {
  const rows: ScatterRow[] = [];
  for (const ix of indices) {
    if (ix.basis <= 0) continue;
    for (const pk of ix.peaks) {
      const pred = ix.predicted_q[pk.ratio_position - 1];
      if (pred == null) continue;
      const ratio = pred / ix.basis;
      rows.push({ ratio, q: pk.q_observed, phase: ix.phase, color: phaseColor(ix.phase) });
    }
  }
  return rows;
}

export function MillerPlot({ indices }: MillerPlotProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const data = toScatterData(indices);

    const el = Plot.plot({
      width:  host.clientWidth  || 360,
      height: host.clientHeight || 260,
      marginLeft: 50, marginBottom: 40,
      x: { label: "√(h²+k²+l²) (normalized)" },
      y: { label: "q (nm⁻¹)" },
      marks: data.length === 0 ? [] : [
        Plot.linearRegressionY(data, {
          x: "ratio", y: "q",
          stroke: "var(--color-fg-muted)",
          strokeDasharray: "4,3",
        }),
        Plot.dot(data, {
          x: "ratio", y: "q",
          fill: (d: ScatterRow) => d.color,
          stroke: "var(--color-bg)",
          strokeWidth: 1,
          r: 4,
          title: (d: ScatterRow) => `${d.phase}\nratio ${d.ratio.toFixed(3)}\nq ${d.q.toFixed(4)}`,
        }),
      ],
    });

    host.replaceChildren(el);
    return () => { host.replaceChildren(); };
  }, [indices]);

  return <div ref={hostRef} className="w-full h-full" data-testid="miller-plot" />;
}
