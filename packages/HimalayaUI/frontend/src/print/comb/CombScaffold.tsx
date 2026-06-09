import { type ReactNode } from "react";
import { makeAxis, axisTicks, type Axis1D } from "../plot/projection";

export const GUTTER_W = 96;
export const ROW_H = 48;
export const TOP_PAD = 10;
export const AXIS_H = 28;
export const PLOT_PAD_X = 14;
export const MIN_PX_PER_DECADE = 240;
export const MIN_PLOT_W = 320;
const DEFAULT_MAX_W = 760;

export interface ScaffoldRow {
  gutterTitle: string;
  gutterSub?: string;
  preview?: boolean;
}

/** What the scaffold hands the render-prop so a chart can place glyphs. */
export interface ScaffoldCtx {
  /** Log-q axis mapping data → pixel across the plot area. */
  x: Axis1D;
  plotW: number;
  padX: number;
  rowH: number;
  topPad: number;
  rowCount: number;
  /** y of row i's baseline (teeth/rings/zero-line sit here; teeth grow upward). */
  baselineY: (i: number) => number;
  rowsBottom: number;
}

interface Props {
  rows: ScaffoldRow[];
  xDomain: [number, number];
  maxWidth?: number;
  ariaLabel: string;
  children: (ctx: ScaffoldCtx) => ReactNode;
}

/** Compact q tick label: 0.05, 0.1, 0.12, 0.2 … */
function formatTick(q: number): string {
  return String(Number(q.toPrecision(2)));
}

export function CombScaffold({ rows, xDomain, maxWidth, ariaLabel, children }: Props): JSX.Element {
  const maxW = maxWidth ?? DEFAULT_MAX_W;
  const [lo, hi] = xDomain;
  // Always log. Floor the plot width so √N teeth never crowd below legibility;
  // when the container is narrower than gutter+floor, the pane scrolls.
  const decades = Math.max(0.5, Math.log10(hi / lo));
  const floorPlotW = Math.max(MIN_PLOT_W, Math.ceil(decades * MIN_PX_PER_DECADE));
  const fitPlotW = Math.max(0, maxW - GUTTER_W - 2 * PLOT_PAD_X);
  const plotW = Math.max(floorPlotW, fitPlotW);
  const svgW = plotW + 2 * PLOT_PAD_X;
  const rowsBottom = TOP_PAD + rows.length * ROW_H;
  const svgH = rowsBottom + AXIS_H;

  const x = makeAxis([lo, hi], [PLOT_PAD_X, PLOT_PAD_X + plotW], "log");
  const baselineY = (i: number): number => TOP_PAD + i * ROW_H + ROW_H - 14;
  const ctx: ScaffoldCtx = {
    x, plotW, padX: PLOT_PAD_X, rowH: ROW_H, topPad: TOP_PAD,
    rowCount: rows.length, baselineY, rowsBottom,
  };
  const ticks = axisTicks(x, 6);

  return (
    <div className="flex items-stretch" style={{ maxWidth: maxW }} data-testid="comb-scaffold">
      {/* Pinned gutter (HTML, never scrolled) — vertically aligned to the SVG rows.
          The top pad is the container's own padding (NOT a margin on row 0): a child
          margin-top would collapse through this block — which has no BFC — and shove
          the whole gutter down relative to the SVG pane (an overflow BFC) whose row
          math already bakes in TOP_PAD. Padding the container keeps the two aligned. */}
      <div className="shrink-0 select-none" style={{ width: GUTTER_W, paddingTop: TOP_PAD }}>
        {rows.map((r, i) => (
          <div
            key={i}
            data-role="gutter-row"
            {...(r.preview ? { "data-preview": "true" } : {})}
            className="flex flex-col justify-end pr-2"
            style={{ height: ROW_H, paddingBottom: 8 }}
          >
            <div className={`font-mono text-data leading-none ${r.preview ? "text-ink-soft" : "text-ink"}`}>
              {r.gutterTitle}
            </div>
            {r.gutterSub ? (
              <div className="font-mono text-meta text-ink-soft leading-none mt-0.5">{r.gutterSub}</div>
            ) : null}
          </div>
        ))}
      </div>

      {/* Scrollable q-pane. */}
      <div className="overflow-x-auto" style={{ maxWidth: maxW - GUTTER_W }}>
        <svg width={svgW} height={svgH} role="img" aria-label={ariaLabel} data-plot-w={plotW}>
          {ticks.map((t, i) => {
            const px = x.to(t.value);
            return (
              <g key={i} data-role="q-tick" data-tick-kind={t.kind}>
                <line
                  x1={px} y1={TOP_PAD} x2={px} y2={rowsBottom}
                  stroke="var(--color-hair)" strokeWidth={1}
                  opacity={t.kind === "major" ? 0.9 : 0.5}
                />
                {t.kind !== "minor" ? (
                  <text
                    x={px} y={rowsBottom + 16} textAnchor="middle"
                    className="font-mono" fontSize={9.5} fill="var(--color-ink-faint)"
                  >
                    {formatTick(t.value)}
                  </text>
                ) : null}
              </g>
            );
          })}
          {rows.map((r, i) => (
            <line
              key={i}
              data-role="row-baseline"
              {...(r.preview ? { "data-preview": "true" } : {})}
              x1={PLOT_PAD_X} y1={baselineY(i)} x2={PLOT_PAD_X + plotW} y2={baselineY(i)}
              stroke="var(--color-hair-strong)" strokeWidth={1}
              {...(r.preview ? { strokeDasharray: "3 3" } : {})}
              opacity={0.8}
            />
          ))}
          {children(ctx)}
        </svg>
      </div>
    </div>
  );
}
