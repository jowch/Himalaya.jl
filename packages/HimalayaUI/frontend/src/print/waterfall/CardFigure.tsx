import { useEffect, useRef, useState } from "react";
import { TracePlot, type TraceModel } from "../plot/TracePlot";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { waterfallQDomain, type WaterfallRow } from "./waterfallModel";
import { phaseColor } from "../../phases";

// Frozen mini-waterfall geometry — mirrors the SeriesCard `#wf-*` figure in
// docs/redesign-mockups/series-folio.html (W=340 reference, stepY=31, pads).
const PAD_L = 14;
const PAD_R = 14;
const PAD_T = 12;
const PAD_B = 12;
const STEP_Y = 31;
const Y_HEADROOM = 0.12; // match WaterfallChart's per-box headroom

/** Mockup reference width — used until the container is measured (and in
 *  non-DOM tests, where ResizeObserver is a no-op stub). */
const DEFAULT_WIDTH = 340;

const UNINDEXED_COLOR = "var(--color-ink-faint)";

export interface CardFigureProps {
  /** Rows low→high (display order); rendered bottom-up. */
  rows: WaterfallRow[];
  /** Fixed outer width in px. Omit (the default) to fill the container and
   *  track its width (FOL-FIGCLIP — a ~301px card must not clip the high-q
   *  tail of a 340px figure). */
  width?: number;
  /** PLACEMENT ONLY. */
  className?: string;
}

function anchorsToPeaks(row: WaterfallRow): PlotPeak[] {
  return row.anchors.map((a) => ({ id: a.id, q: a.q, intensity: a.intensity, source: "auto" as const }));
}

/**
 * A frozen, non-interactive miniature of {@link WaterfallChart} for embedding in
 * a series card. No axis, no label gutter, no hover — just stacked phase-coloured
 * trace rows over a shared log-q domain. Row 0 sits at the bottom.
 *
 * Width-responsive by default: with no `width` prop the figure fills its
 * container and re-renders the rows at the measured width (same
 * ResizeObserver pattern as PlotFrame's controlled-width fallback).
 */
export function CardFigure({ rows, width, className = "" }: CardFigureProps): JSX.Element {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [measured, setMeasured] = useState<number | null>(null);

  // Measure the container unless an explicit width pins the figure.
  useEffect(() => {
    if (width !== undefined) return;
    const el = containerRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width ?? 0;
      if (w > 0) setMeasured(w);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [width]);

  const w = width ?? measured ?? DEFAULT_WIDTH;
  const qDomain = waterfallQDomain(rows);
  const pw = Math.max(0, w - PAD_L - PAD_R);
  const height = PAD_T + rows.length * STEP_Y + PAD_B;

  return (
    <div
      ref={containerRef}
      className={`relative ${className}`}
      style={{ width: width ?? "100%", height }}
      data-testid="card-figure"
      data-row-count={rows.length}
    >
      {rows.map((row, i) => {
        const model: TraceModel = { trace: row.trace, peaks: anchorsToPeaks(row), phase: row.phase };
        // Bottom-up: row 0 at the bottom. Top of band i counts up from the base.
        const top = PAD_T + (rows.length - 1 - i) * STEP_Y;
        const color = row.phase ? phaseColor(row.phase) : UNINDEXED_COLOR;
        return (
          <div
            key={row.key}
            data-row-key={row.key}
            data-phase={row.phase ?? "none"}
            className="absolute"
            style={{ left: PAD_L, top, width: pw, height: STEP_Y, color }}
          >
            <TracePlot
              trace={model}
              axes={false}
              interaction={false}
              xType="log"
              xDomain={qDomain}
              yHeadroom={Y_HEADROOM}
              height={STEP_Y}
              width={pw}
              paperColor="var(--color-plate)"
              show={{ peaks: true, labels: false, band: false }}
            />
          </div>
        );
      })}
    </div>
  );
}
