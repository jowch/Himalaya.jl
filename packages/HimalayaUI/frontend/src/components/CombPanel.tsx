import { useState } from "react";
import { useAppState } from "../state";
import { usePeaks, useIndices, useAssignment } from "../queries";
import { deriveActiveIndices } from "../lib/assignment";
import { phaseColor } from "../phases";
import { HintText } from "./ui";
import { SegmentedControl } from "./ui/SegmentedControl";
import type { IndexEntry, Peak } from "../api";

/**
 * CombPanel (Plan D Task D-5) — predicted reflections as comb teeth on a shared
 * log-q ruler; the q-link triple's third surface.
 *
 *  - One row of teeth per assigned phase, drawn from `IndexEntry.predicted_q`
 *    (all orders). A predicted order that matches an observed peak is a solid
 *    tooth; a predicted-but-absent order is a hollow caret.
 *  - A leftover row collects observed peaks that no assigned phase claims.
 *  - A SegmentedControl toggles to the indexing-space residual view
 *    (Δq/q vs predicted q; a correct phase sits on the zero line).
 *  - A dashed ghost row previews the hovered candidate (Plan D-7).
 *
 * q-link (Plan C §5.1, LOW-2): the tooth whose q matches the global `hoveredQ`
 * is HOT — it grows and gains a terracotta ring (data-hot="true"), NOT a
 * recolour. Hovering a tooth pushes its q back through `setHoveredQ`, so the
 * trace peak / detector ring / comb tooth light together.
 *
 * Visual contract: docs/redesign-mockups/2026-05-29-focus-plot.html (#combs).
 * Drawn as a dedicated SVG ruler; phase colour threads through SVG attributes
 * (design-guard-safe), className is placement-only.
 */

const VB = { l: 96, r: 22, t: 34, b: 30, w: 620, h: 300 };
const PLOT_W = VB.w - VB.l - VB.r;
const ACCENT = "var(--color-accent)";

function logX(q: number, qMin: number, qMax: number): number {
  if (q <= 0 || qMin <= 0 || qMax <= qMin) return VB.l;
  return VB.l + ((Math.log(q) - Math.log(qMin)) / (Math.log(qMax) - Math.log(qMin))) * PLOT_W;
}

const COMB_VIEWS = [
  { value: "comb" as const, label: "comb" },
  { value: "resid" as const, label: "residual" },
];

export function CombPanel(): JSX.Element {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const hoveredQ = useAppState((s) => s.hoveredQ);
  const setHoveredQ = useAppState((s) => s.setHoveredQ);
  const previewIndexId = useAppState((s) => s.previewIndexId);
  const [view, setView] = useState<"comb" | "resid">("comb");

  const peaksQ = usePeaks(activeExposureId);
  const indicesQ = useIndices(activeExposureId);
  const assignmentQ = useAssignment(activeExposureId);

  const peaks: Peak[] = (peaksQ.data ?? []).filter((p) => !p.excluded);
  const indices = indicesQ.data ?? [];
  const active = deriveActiveIndices(assignmentQ.data, indices);
  const previewIndex = previewIndexId != null
    ? indices.find((ix) => ix.id === previewIndexId)
    : undefined;

  if (activeExposureId === undefined) {
    return (
      <div className="flex-1 flex items-center justify-center">
        <HintText>Pick a sample to see its reflections.</HintText>
      </div>
    );
  }

  // log-q domain from the union of observed + predicted q (so absent orders fit).
  const qs = [
    ...peaks.map((p) => p.q),
    ...active.flatMap((ix) => ix.predicted_q),
  ].filter((q) => q > 0);
  const qMin = qs.length ? Math.min(...qs) * 0.9 : 0.01;
  const qMax = qs.length ? Math.max(...qs) * 1.1 : 0.2;

  // q-link tolerance — span-relative (mirror the ring rule).
  const peakQs = peaks.map((p) => p.q);
  const qLo = peakQs.length ? Math.min(...peakQs) : 0;
  const qHi = peakQs.length ? Math.max(...peakQs) : 1;
  const tol = Math.max((qHi - qLo) * 0.02, 1e-6);
  const isHot = (q: number): boolean =>
    hoveredQ !== undefined && Math.abs(q - hoveredQ) <= tol;

  // Which observed peak (if any) sits at a predicted q — claim by tolerance.
  const claimedPeakIds = new Set<number>();
  for (const ix of active) for (const p of ix.peaks) claimedPeakIds.add(p.peak_id);
  const leftover = peaks.filter((p) => !claimedPeakIds.has(p.id));
  const nearestPeakQ = (q: number): number | undefined => {
    let best: Peak | undefined; let bestD = Infinity;
    for (const p of peaks) {
      const d = Math.abs(p.q - q);
      if (d < bestD) { bestD = d; best = p; }
    }
    return best && bestD <= tol ? best.q : undefined;
  };

  const railY = VB.h - VB.b;
  const top = VB.t;

  const tools = (
    <SegmentedControl
      options={COMB_VIEWS}
      value={view}
      onChange={setView}
      role="radiogroup"
      size="sm"
      aria-label="Comb view"
      testId="comb-view"
      className="self-end"
    />
  );

  // ── residual view ──
  if (view === "resid") {
    const midY = (top + railY) / 2;
    const span = (railY - top) / 2 - 8;
    const yResid = (d: number): number => midY - Math.max(-1, Math.min(1, d / 0.03)) * span;
    return (
      <div className="flex flex-col h-full min-h-0">
        <div className="flex items-center justify-between px-1 pb-1">
          <span className="text-xs font-semibold text-ink-soft">Reflections — residual</span>
          {tools}
        </div>
        <div className="flex-1 min-h-0">
          <svg data-testid="comb-svg" viewBox={`0 0 ${VB.w} ${VB.h}`}
               preserveAspectRatio="xMidYMid meet" className="block w-full h-full"
               aria-label="Indexing-space residual: Δq/q vs predicted q">
            <line data-testid="comb-residual-zero-line"
                  x1={VB.l} y1={midY} x2={VB.w - VB.r} y2={midY}
                  stroke="var(--color-hair-strong)" strokeWidth={1} />
            <text x={8} y={midY + 3} fill="var(--color-ink-faint)" fontSize={9.5}>Δq/q = 0</text>
            <text x={8} y={top + 8} fill="var(--color-ink-faint)" fontSize={9.5}>residual</text>
            {active.flatMap((ix) =>
              ix.peaks.map((p) => {
                if (p.q_observed <= 0) return null;
                // Δq/q against predicted-q (= observed − residual). Guard a
                // near-zero predicted-q (a degenerate index) explicitly so dd
                // stays bounded rather than relying solely on yResid's clamp.
                const predictedQ = p.q_observed - p.residual;
                const dd = p.residual / (Math.abs(predictedQ) > 1e-9 ? predictedQ : p.q_observed);
                const x = logX(p.q_observed, qMin, qMax);
                const y = yResid(dd);
                const color = phaseColor(ix.phase);
                const hot = isHot(p.q_observed);
                return (
                  <g key={`r-${ix.id}-${p.ratio_position}`}
                     data-testid={`comb-resid-${ix.id}-${p.ratio_position}`}
                     data-hot={hot || undefined}>
                    <line x1={x} y1={midY} x2={x} y2={y}
                          stroke={hot ? ACCENT : color} strokeWidth={hot ? 2 : 1} opacity={0.6} />
                    <circle cx={x} cy={y} r={hot ? 4 : 2.6}
                            fill={hot ? "none" : color} stroke={hot ? ACCENT : "none"}
                            strokeWidth={hot ? 2 : 0} />
                  </g>
                );
              }))}
          </svg>
        </div>
      </div>
    );
  }

  // ── comb view ──
  const combRows: { ix: IndexEntry; ghost: boolean }[] = active.map((ix) => ({ ix, ghost: false }));
  if (previewIndex && !active.some((a) => a.id === previewIndex.id)) {
    combRows.push({ ix: previewIndex, ghost: true });
  }
  const nRows = combRows.length + 1; // + leftover row
  const band = railY - top;
  const rowGap = band / nRows;
  const toothH = Math.min(rowGap * 0.58, 42);

  return (
    <div className="flex flex-col h-full min-h-0">
      <div className="flex items-center justify-between px-1 pb-1">
        <span className="text-xs font-semibold text-ink-soft">Reflections — comb</span>
        {tools}
      </div>
      <div className="flex-1 min-h-0">
        <svg data-testid="comb-svg" viewBox={`0 0 ${VB.w} ${VB.h}`}
             preserveAspectRatio="xMidYMid meet" className="block w-full h-full"
             aria-label="Predicted reflections as comb teeth on a shared q ruler">
          {/* shared q ruler */}
          <line x1={VB.l} y1={railY} x2={VB.w - VB.r} y2={railY}
                stroke="var(--color-hair-strong)" strokeWidth={1} />
          {/* observed-peak vertical guides */}
          {peaks.map((p) => {
            const x = logX(p.q, qMin, qMax);
            return (
              <line key={`g-${p.id}`} x1={x} y1={top} x2={x} y2={railY}
                    stroke="var(--color-hair)" strokeWidth={1}
                    strokeDasharray="1 3" opacity={0.8} />
            );
          })}

          {combRows.map(({ ix, ghost }, ri) => {
            const yBase = top + (ri + 0.72) * rowGap;
            const color = phaseColor(ix.phase);
            const claimedQ = new Set(ix.peaks.map((p) => p.q_observed));
            return (
              <g key={`row-${ix.id}`}>
                <text className="font-bold" x={8} y={yBase - 8} fill={color}
                      fontSize={11} opacity={ghost ? 0.8 : 1}>
                  {ghost ? `+ ${ix.phase}` : ix.phase}
                </text>
                <text x={8} y={yBase + 6} fill="var(--color-ink-faint)" fontSize={9.5}>
                  {ghost ? "hypothetical" : (ix.lattice_d != null ? `a = ${ix.lattice_d.toFixed(0)}` : "")}
                </text>
                <line x1={VB.l} y1={yBase} x2={VB.w - VB.r} y2={yBase}
                      stroke={ghost ? "var(--color-hair)" : "var(--color-hair-strong)"}
                      strokeWidth={1} strokeDasharray={ghost ? "3 3" : undefined} />
                {ix.predicted_q.map((q, oi) => {
                  if (q < qMin || q > qMax) return null;
                  const order = oi + 1;
                  const x = logX(q, qMin, qMax);
                  // predicted order is "observed" if it claims a peak (a member
                  // q_observed within tolerance) — else it is absent.
                  const matchQ = nearestPeakQ(q);
                  const observed = matchQ !== undefined
                    && [...claimedQ].some((cq) => Math.abs(cq - matchQ) <= tol);
                  const hot = !ghost && matchQ !== undefined && isHot(matchQ);
                  if (observed) {
                    return (
                      <g key={`t-${ix.id}-${order}`}
                         data-testid={`comb-tooth-${ix.id}-${order}`}
                         data-hot={hot || undefined}
                         style={{ cursor: "pointer" }}
                         onMouseEnter={() => !ghost && matchQ !== undefined && setHoveredQ(matchQ)}
                         onMouseLeave={() => !ghost && setHoveredQ(undefined)}>
                        <line x1={x} y1={yBase} x2={x} y2={yBase - (hot ? toothH * 1.12 : toothH)}
                              stroke={hot ? ACCENT : color}
                              strokeWidth={hot ? 2.8 : ghost ? 1.6 : 2}
                              strokeDasharray={ghost ? "2.5 2" : undefined} />
                        <circle cx={x} cy={yBase - (hot ? toothH * 1.12 : toothH)}
                                r={hot ? 5 : 2.4}
                                fill={ghost ? "var(--color-plate)" : hot ? "none" : color}
                                stroke={hot ? ACCENT : ghost ? color : "none"}
                                strokeWidth={hot ? 2 : ghost ? 1.4 : 0} />
                        <text x={x} y={yBase - toothH - 6} textAnchor="middle"
                              fontSize={9.5} fontWeight={700}
                              fill={hot ? ACCENT : color}
                              className="font-mono">
                          {String(order)}
                        </text>
                        {/* generous transparent hit target */}
                        <circle cx={x} cy={yBase - toothH / 2} r={10} fill="transparent" />
                      </g>
                    );
                  }
                  // predicted, absent → hollow caret
                  const ah = toothH * 0.72;
                  return (
                    <g key={`a-${ix.id}-${order}`}
                       data-testid={`comb-absent-${ix.id}-${order}`}>
                      <line x1={x} y1={yBase} x2={x} y2={yBase - ah}
                            stroke="var(--color-ink-faint)" strokeWidth={1.3}
                            strokeDasharray="1.5 1.8" />
                      <path d={`M${x - 3.2} ${yBase - ah} L${x} ${yBase - ah - 4} L${x + 3.2} ${yBase - ah}`}
                            fill="none" stroke="var(--color-ink-faint)" strokeWidth={1.3} />
                    </g>
                  );
                })}
              </g>
            );
          })}

          {/* leftover row */}
          {(() => {
            const yL = top + (combRows.length + 0.55) * rowGap;
            return (
              <g>
                <text x={8} y={yL + 3} fill="var(--color-ink-faint)" fontSize={9.5}>
                  {leftover.length ? "leftover" : "leftover — none"}
                </text>
                {leftover.map((p) => {
                  const x = logX(p.q, qMin, qMax);
                  const hot = isHot(p.q);
                  return (
                    <g key={`l-${p.id}`} data-testid={`comb-leftover-${p.id}`}
                       data-hot={hot || undefined}
                       style={{ cursor: "pointer" }}
                       onMouseEnter={() => setHoveredQ(p.q)}
                       onMouseLeave={() => setHoveredQ(undefined)}>
                      <circle cx={x} cy={yL} r={hot ? 5 : 3.6} fill="none"
                              stroke={hot ? ACCENT : "var(--color-ink-faint)"}
                              strokeWidth={hot ? 2.4 : 1.6} />
                      <circle cx={x} cy={yL} r={10} fill="transparent" />
                    </g>
                  );
                })}
              </g>
            );
          })()}
        </svg>
      </div>
    </div>
  );
}
