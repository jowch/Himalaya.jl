import { type ReactNode, useLayoutEffect, useRef } from "react";
import { makeAxis, axisTicks, type Axis1D } from "../plot/projection";
import { formatAxis } from "../../lib/plot/formatAxis";

/** Low-q lead-in (px) kept visible before the first reflection when the comb is
 *  scrolled to the peaks, so the first tooth isn't flush against the gutter. */
const FOCUS_LEAD_IN = 40;

export const GUTTER_W = 96;
export const ROW_H = 48;
export const TOP_PAD = 10;
export const AXIS_H = 44; // tick numbers + the "q (Å⁻¹)" axis title below them
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
  /** q-extent of the reflections. When the pane overflows, it default-scrolls to
   *  show these (centred if they fit, else the first reflection near the left)
   *  instead of the low-q beam dropoff. Runs on mount + when this/the domain/width
   *  change — never on hover or after a manual scroll (none change these deps). */
  focusQRange?: [number, number];
  maxWidth?: number;
  ariaLabel: string;
  children: (ctx: ScaffoldCtx) => ReactNode;
}

export function CombScaffold({ rows, xDomain, focusQRange, maxWidth, ariaLabel, children }: Props): JSX.Element {
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

  // FO-COMBSCROLL-PEAKS: default-scroll the q-pane to the reflections so the panel
  // opens on the peaks, not the low-q beam dropoff (which reads as empty until the
  // user discovers they can scroll). Deps are PRIMITIVES (the focus q's, the plot
  // width, the domain) so this fires on mount and when the content/geometry change
  // — never on hover (no dep changes) or after a manual scroll (scrolling fires no
  // render). x.to is pure given lo/hi/plotW, so omitting `x` from deps is safe.
  const paneRef = useRef<HTMLDivElement>(null);
  const focusLo = focusQRange?.[0];
  const focusHi = focusQRange?.[1];
  useLayoutEffect(() => {
    const pane = paneRef.current;
    if (!pane || focusLo === undefined || focusHi === undefined) return;
    const xFirst = x.to(focusLo);
    const xLast = x.to(focusHi);
    if (!Number.isFinite(xFirst) || !Number.isFinite(xLast)) return;
    const paneW = pane.clientWidth;
    const span = xLast - xFirst;
    // Centre the reflections when they fit; otherwise put the first one near the
    // left so the comb reads from its start.
    const target =
      span <= paneW - 2 * FOCUS_LEAD_IN
        ? (xFirst + xLast) / 2 - paneW / 2
        : xFirst - FOCUS_LEAD_IN;
    pane.scrollLeft = Math.max(0, target);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [focusLo, focusHi, plotW, lo, hi]);

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

      {/* Scrollable q-pane. FO-COMB-SCROLL: when the q-axis is wider than the
          pane it scrolls horizontally; label it and make it keyboard-focusable
          so the affordance is not a cryptic unlabelled scrollbar. */}
      <div
        ref={paneRef}
        className="overflow-x-auto"
        style={{ maxWidth: maxW - GUTTER_W }}
        tabIndex={0}
        role="group"
        aria-label="Reflection comb q-axis, scroll horizontally"
      >
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
                  // FO-COMB-AXIS / CC-AXIS-TICK: the q-values are the data, so the
                  // tick numbers read at ink-soft (AA), not the faint axis tone.
                  <text
                    x={px} y={rowsBottom + 16} textAnchor="middle"
                    className="font-mono" fontSize={9.5} fill="var(--color-ink-soft)"
                  >
                    {formatAxis(t.value)}
                  </text>
                ) : null}
              </g>
            );
          })}
          {/* FO-COMB-AXIS: name the axis so the comb reads as a quantitative
              figure (it shares the trace's q-domain, so a q here is the same q
              there). Centred under the plot area. */}
          <text
            data-role="comb-axis-title"
            x={PLOT_PAD_X + plotW / 2}
            y={rowsBottom + 36}
            textAnchor="middle"
            className="font-sans"
            fontSize={11}
            fontWeight={600}
            fill="var(--color-ink-soft)"
          >
            q (Å⁻¹)
          </text>
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
