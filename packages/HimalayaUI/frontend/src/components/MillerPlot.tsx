import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { IndexEntry } from "../api";
import { phaseColor } from "../phases";

export interface MillerPlotProps {
  indices: IndexEntry[];
  /** When set, dims all indices except this one. Also accepts candidates not in `indices`. */
  hoveredIndex?: IndexEntry | undefined;
}

export interface ScatterRow {
  ratio: number;
  q: number;
  phase: string;
  color: string;
  /** Stable grouping key so the regression fits one line per index, not one for all. */
  indexId: number;
}

export function toScatterData(indices: IndexEntry[]): ScatterRow[] {
  const rows: ScatterRow[] = [];
  for (const ix of indices) {
    if (ix.basis <= 0) continue;
    for (const pk of ix.peaks) {
      const pred = ix.predicted_q[pk.ratio_position - 1];
      if (pred == null) continue;
      const ratio = pred / ix.basis;
      rows.push({
        ratio,
        q: pk.q_observed,
        phase: ix.phase,
        color: phaseColor(ix.phase),
        indexId: ix.id,
      });
    }
  }
  return rows;
}

export function MillerPlot({ indices, hoveredIndex }: MillerPlotProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const hoveredId = hoveredIndex?.id;
    const dimming = hoveredId !== undefined;

    // Combine active indices with the hovered candidate (if it's not already active).
    const activeIds = new Set(indices.map((ix) => ix.id));
    const allIndices = hoveredIndex && !activeIds.has(hoveredIndex.id)
      ? [...indices, hoveredIndex]
      : indices;

    const data = toScatterData(allIndices);

    const byIndex = new Map<number, ScatterRow[]>();
    for (const r of data) {
      const list = byIndex.get(r.indexId);
      if (list) list.push(r); else byIndex.set(r.indexId, [r]);
    }

    const regressionMarks: Plot.Markish[] = [];
    for (const [id, rows] of byIndex) {
      if (rows.length < 2) continue;
      const color = rows[0]!.color;
      const isHovered = id === hoveredId;
      const opacity = dimming ? (isHovered ? 1 : 0.18) : 0.85;
      regressionMarks.push(
        Plot.linearRegressionY(rows, {
          x: "ratio",
          y: "q",
          stroke: color,
          strokeOpacity: opacity,
          strokeWidth: isHovered ? 1.5 : 1,
          ...(isHovered ? {} : { strokeDasharray: "4,3" }),
        }),
      );
    }

    const el = Plot.plot({
      width:  host.clientWidth  || 360,
      height: host.clientHeight || 260,
      marginLeft: 32, marginBottom: 22, marginTop: 6, marginRight: 8,
      style: {
        fontFamily: "var(--font-sans)",
        color: "var(--color-fg-muted)",
        background: "transparent",
        overflow: "visible",
        fontSize: "9px",
      },
      x: { label: null, ticks: 4 },
      y: { label: null, ticks: 3 },
      marks: data.length === 0 ? [] : [
        ...regressionMarks,
        Plot.dot(data, {
          x: "ratio", y: "q",
          fill: (d: ScatterRow) => d.color,
          fillOpacity: dimming
            ? (d: ScatterRow) => (d.indexId === hoveredId ? 1 : 0.15)
            : 1,
          stroke: "var(--color-bg)",
          strokeWidth: 1,
          r: dimming
            ? (d: ScatterRow) => (d.indexId === hoveredId ? 4 : 2.5)
            : 3,
          title: (d: ScatterRow) => `${d.phase}\nratio ${d.ratio.toFixed(3)}\nq ${d.q.toFixed(4)}`,
        }),
      ],
    });

    host.replaceChildren(el);
    return () => { host.replaceChildren(); };
  }, [indices, hoveredIndex]);

  return <div ref={hostRef} className="w-full h-full" data-testid="miller-plot" />;
}
