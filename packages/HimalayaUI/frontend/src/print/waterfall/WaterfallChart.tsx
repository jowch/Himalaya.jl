import { useEffect, useRef, useState } from "react";
import { TracePlot, type TraceModel } from "../plot/TracePlot";
import { Axis } from "../plot/Axis";
import { makeAxis } from "../plot/projection";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { waterfallQDomain, type WaterfallRow } from "./waterfallModel";
import { snapToPeakQ } from "./cursor";
import { PEAK_HIT_PX } from "../plot/interaction";

const AXIS_H = 44;          // height of the shared bottom axis strip
const LABEL_W = 56;         // right-margin label gutter
const TOTAL_H = 420;        // fixed internal stack height (CSS scales width)
const Y_HEADROOM = 0.12;    // each box keeps ~12% headroom above its peak

export interface WaterfallChartProps {
  /** Rows low→high (display order); rendered bottom-up. */
  rows: WaterfallRow[];
  xType?: "log" | "linear";
  /** Controlled hot row; falls back to internal hover when omitted. */
  hoveredKey?: string;
  onHoverRow?: (key?: string) => void;
  onHoverQ?: (q?: number) => void;
  /** Controlled cursor q; falls back to internal hover when omitted. */
  hoveredQ?: number;
  /** Fit-to-width ceiling (px). Default 1080 (mockup plate). */
  maxWidth?: number;
  /** CSS min-width floor on the figure; it won't shrink below this. Default 560. */
  minWidth?: number;
  /** PLACEMENT ONLY. */
  className?: string;
}

function anchorsToPeaks(row: WaterfallRow): PlotPeak[] {
  return row.anchors.map((a) => ({ id: a.id, q: a.q, intensity: a.intensity, source: "auto" as const }));
}

export function WaterfallChart({
  rows,
  xType = "log",
  hoveredKey,
  onHoverRow,
  onHoverQ,
  hoveredQ,
  maxWidth = 1080,
  minWidth = 560,
  className = "",
}: WaterfallChartProps): JSX.Element {
  const [internalHot, setInternalHot] = useState<string | undefined>(undefined);
  const hot = hoveredKey ?? internalHot;
  const [internalHoveredQ, setInternalHoveredQ] = useState<number | undefined>(undefined);
  const cursorQ = hoveredQ ?? internalHoveredQ;

  const rootRef = useRef<HTMLDivElement | null>(null);
  const [measuredW, setMeasuredW] = useState<number | null>(null);

  useEffect(() => {
    const el = rootRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width ?? 0;
      if (w > 0) setMeasuredW(w);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const qDomain = waterfallQDomain(rows);

  // Effective figure width: the measured container (CSS-clamped to [minWidth, maxWidth]),
  // falling back to maxWidth before first measure / in non-DOM tests.
  const effectiveW = Math.min(measuredW ?? maxWidth, maxWidth);
  const plotW = Math.max(0, effectiveW - LABEL_W);

  const sharedX = makeAxis(qDomain, [4, plotW - 4], xType);

  const totalWeight = rows.reduce((s, r) => s + Math.max(0, r.bandHeight), 0) || rows.length;

  // Bottom-up: display-order-0 (rows[0]) sits at the bottom. Walk rows in
  // reverse for top→bottom pixel placement.
  let cumulative = 0;
  const placed = [...rows].reverse().map((row) => {
    const h = (Math.max(0, row.bandHeight) / totalWeight) * TOTAL_H;
    // TODO: clamp/scale yOffset into TOTAL_H so a large durable offset can't push a row off-canvas.
    const top = cumulative + row.yOffset;
    cumulative += h;
    return { row, top, h };
  });

  const setHot = (key?: string): void => {
    if (hoveredKey === undefined) setInternalHot(key);
    onHoverRow?.(key);
  };

  const setCursorQ = (q: number | undefined): void => {
    if (hoveredQ === undefined) setInternalHoveredQ(q);
    onHoverQ?.(q);
  };

  return (
    <div ref={rootRef} className={`relative w-full ${className}`} style={{ maxWidth, minWidth }} data-testid="waterfall" data-xtype={xType}>
      <div
        className="relative"
        style={{ height: TOTAL_H }}
        data-testid="wf-stack"
        onPointerMove={(e) => {
          const rect = e.currentTarget.getBoundingClientRect();
          const snapped = snapToPeakQ(e.clientX - rect.left, rows, (q) => sharedX.to(q), PEAK_HIT_PX);
          setCursorQ(snapped ?? undefined);
        }}
        onPointerLeave={() => setCursorQ(undefined)}
      >
        {placed.map(({ row, top, h }) => {
          const model: TraceModel = { trace: row.trace, peaks: anchorsToPeaks(row), phase: row.phase };
          const isHot = hot === row.key;
          const isDim = hot !== undefined && !isHot;
          return (
            <div
              key={row.key}
              data-role="wf-row"
              data-key={row.key}
              {...(isHot ? { "data-hot": "true" } : {})}
              {...(isDim ? { "data-dim": "true" } : {})}
              className="absolute left-0 flex items-center"
              style={{ top, height: h, width: "100%", opacity: isDim ? 0.32 : 1 }}
              onMouseEnter={() => setHot(row.key)}
              onMouseLeave={() => setHot(undefined)}
            >
              <TracePlot
                trace={model}
                axes={false}
                xType={xType}
                xDomain={qDomain}
                yHeadroom={Y_HEADROOM}
                height={h}
                width={plotW}
                paperColor="var(--color-plate)"
                show={{ peaks: true, labels: false, band: false }}
              />
              <div
                data-role="wf-label"
                className={`text-meta font-mono leading-none pl-2 ${isHot ? "text-ink" : "text-ink-soft"}`}
                style={{ width: LABEL_W }}
              >
                {row.label}
              </div>
            </div>
          );
        })}
        {cursorQ != null ? (
          <>
            <div
              data-role="wf-qguide"
              data-testid="wf-qguide"
              data-q={cursorQ}
              className="absolute top-0 w-px bg-accent pointer-events-none"
              style={{ left: sharedX.to(cursorQ), height: TOTAL_H }}
            />
            <div
              data-role="wf-qreadout"
              data-testid="wf-qreadout"
              // Anchored at the FOOT of the guide (bottom edge at the stack base, just
              // above the q-axis) — parallels TracePlot's q-readout sitting at the axis
              // baseline, so the value reads where the eye reads q.
              className="absolute -translate-x-1/2 -translate-y-full bg-plate border border-hair rounded-sm px-1 text-meta font-mono text-ink pointer-events-none"
              style={{ left: sharedX.to(cursorQ), top: TOTAL_H }}
            >
              {cursorQ.toFixed(3)}
            </div>
          </>
        ) : null}
      </div>

      <svg width={plotW} height={AXIS_H} role="img" aria-label="scattering vector q" data-testid="wf-axis">
        {/* Axis baseline sits at the strip top (plotHeight=0); ticks + label descend into the 44px strip. */}
        {/* range [4, plotW-4] matches a row TracePlot's axes=false 4px inset, so ticks sit under the traces. */}
        <Axis axis={sharedX} orientation="bottom" plotWidth={plotW} plotHeight={0} label="q (Å⁻¹)" />
      </svg>
    </div>
  );
}
