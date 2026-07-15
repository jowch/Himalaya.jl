import { useEffect, useRef, useState } from "react";
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
  speculative: boolean;
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
        speculative: ix.kind === "speculative",
      });
    }
  }
  return rows;
}

export function MillerPlot({ indices, hoveredIndex }: MillerPlotProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  const [_resizeKey, setResizeKey] = useState(0);
  useEffect(() => {
    const el = hostRef.current;
    if (!el) return;
    const obs = new ResizeObserver(() => setResizeKey((k) => k + 1));
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

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
      const isHovered    = id === hoveredId;
      const isSpeculative = rows[0]!.speculative;
      const opacity = dimming ? (isHovered ? 1 : 0.18) : 0.85;
      // Speculative regressions are ALWAYS dashed (the visual warning that
      // they're hand-built, sub-minpeaks). For auto, dashes only show when
      // un-hovered as a depth cue.
      const dashed = isSpeculative || !isHovered;
      regressionMarks.push(
        Plot.linearRegressionY(rows, {
          x: "ratio",
          y: "q",
          stroke: color,
          strokeOpacity: opacity,
          strokeWidth: isHovered ? 1.5 : 1,
          // n=2 has zero residual DOF: the confidence band divides by (n-2)
          // and renders an all-NaN <path> (SVG console error on every draw).
          // ci: 0 hides only the band; the regression line always renders.
          ...(rows.length < 3 ? { ci: 0 } : {}),
          ...(dashed ? { strokeDasharray: isSpeculative ? "2,3" : "4,3" } : {}),
        }),
      );
    }

    // Domain is computed from the dot data only — Plot's default inference
    // mixes in `linearRegressionY`'s confidence band, which extrapolates past
    // the data and crops some points off-canvas. Pad 5% on each side so peaks
    // at the extrema sit inside the frame instead of on the axis.
    const xs = data.map((d) => d.ratio);
    const ys = data.map((d) => d.q);
    const padded = (lo: number, hi: number): [number, number] => {
      if (lo === hi) {
        const eps = lo === 0 ? 1 : Math.abs(lo) * 0.05;
        return [lo - eps, hi + eps];
      }
      const pad = (hi - lo) * 0.05;
      return [lo - pad, hi + pad];
    };
    const xDomain = xs.length > 0 ? padded(Math.min(...xs), Math.max(...xs)) : undefined;
    const yDomain = ys.length > 0 ? padded(Math.min(...ys), Math.max(...ys)) : undefined;

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
      x: { label: null, ticks: 4, ...(xDomain ? { domain: xDomain } : {}) },
      y: { label: null, ticks: 3, ...(yDomain ? { domain: yDomain } : {}) },
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
  }, [indices, hoveredIndex, _resizeKey]);

  return <div ref={hostRef} className="w-full h-full" data-testid="miller-plot" />;
}
