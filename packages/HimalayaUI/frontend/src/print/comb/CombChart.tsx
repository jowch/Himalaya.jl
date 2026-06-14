import { CombScaffold, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { dodgeX } from "../plot/marks/PlotLabels";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
const LEFTOVER_COLOR = "var(--color-ink-faint)";
/** Min horizontal gap (px) between adjacent reflection-ratio labels before they
 *  dodge apart (with a leader line back to the tooth). */
const LABEL_MIN_GAP = 20;

interface Props {
  assigned: CombSeries[];
  hovered?: CombSeries;
  leftover: number[];
  /** Incoming q-link; the matching glyph across all rows lights hot. */
  hoveredQ?: number;
  /** Fired on glyph enter (q) / leave (undefined). Omit → inert glyphs. */
  onHoverQ?: (q?: number) => void;
  /** FO-COMB-AXIS: the trace's q-domain. When given, the comb spans the same
   *  q-range as the trace above it (proper multi-tick axis + true reflection
   *  positions) instead of a domain derived only from the reflections present.
   *  Falls back to `combQDomain(rows)` for standalone/story use. */
  xDomain?: [number, number];
  maxWidth?: number;
  className?: string;
}

export function CombChart({ assigned, hovered, leftover, hoveredQ, onHoverQ, xDomain: xDomainProp, maxWidth, className }: Props): JSX.Element {
  const rows = assembleRows(assigned, hovered, leftover);
  const xDomain = xDomainProp ?? combQDomain(rows);
  const focusQRange = reflectionQExtent(rows);
  const scaffoldRows: ScaffoldRow[] = rows.map(rowToGutter);
  return (
    <div className={className}>
      <CombScaffold
        rows={scaffoldRows}
        xDomain={xDomain}
        ariaLabel="reflection comb"
        {...(focusQRange ? { focusQRange } : {})}
        {...(maxWidth !== undefined ? { maxWidth } : {})}
      >
        {(ctx) => rows.map((row, i) => renderRow(row, i, ctx, hoveredQ, onHoverQ))}
      </CombScaffold>
    </div>
  );
}

/** The q-extent of the OBSERVED reflections (solid teeth + leftover rings), so
 *  the scaffold can default-scroll the pane to where the peaks actually are
 *  instead of the low-q beam dropoff. Falls back to ALL teeth when a row has only
 *  predicted-absent carets; `undefined` when there is nothing to scroll to. */
export function reflectionQExtent(rows: CombRow[]): [number, number] | undefined {
  const observed: number[] = [];
  const all: number[] = [];
  for (const row of rows) {
    if (row.kind === "leftover") {
      observed.push(...row.qs);
      all.push(...row.qs);
    } else {
      for (const t of row.series.teeth) {
        all.push(t.q);
        if (t.observed) observed.push(t.q);
      }
    }
  }
  const qs = observed.length > 0 ? observed : all;
  if (qs.length === 0) return undefined;
  return [Math.min(...qs), Math.max(...qs)];
}

function rowToGutter(row: CombRow): ScaffoldRow {
  if (row.kind === "leftover") return { gutterTitle: "leftover", gutterSub: `${row.qs.length} unindexed` };
  const base: ScaffoldRow = {
    gutterTitle: row.series.phase,
    ...(row.series.latticeLabel ? { gutterSub: row.series.latticeLabel } : {}),
  };
  return row.kind === "preview" ? { ...base, preview: true } : base;
}

function renderRow(
  row: CombRow,
  i: number,
  ctx: ScaffoldCtx,
  hoveredQ: number | undefined,
  onHoverQ?: (q?: number) => void,
): JSX.Element {
  const y0 = ctx.baselineY(i);
  const toothH = ctx.rowH - 26;
  const hov = (q: number): { onMouseEnter: () => void; onMouseLeave: () => void } | Record<string, never> =>
    onHoverQ ? { onMouseEnter: () => onHoverQ(q), onMouseLeave: () => onHoverQ(undefined) } : {};
  const isHot = (q: number): boolean => hoveredQ !== undefined && Math.abs(q - hoveredQ) <= TOL;

  if (row.kind === "leftover") {
    return (
      <g key={i} data-role="comb-row" data-row-kind="leftover">
        {row.qs.map((q, j) => {
          const hot = isHot(q);
          return (
            <g
              key={j}
              data-role="leftover-ring"
              data-q={q}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...hov(q)}
            >
              <circle
                cx={ctx.x.to(q)} cy={y0} r={hot ? 5 : 3.6}
                fill="none" stroke={LEFTOVER_COLOR} strokeWidth={hot ? 2.4 : 1.6}
              />
            </g>
          );
        })}
      </g>
    );
  }

  const { color } = row.series;
  // FO-COMBLABEL-DODGE: spread the observed reflection-ratio labels so a dense
  // comb doesn't render them on top of each other; a leader line ties each
  // shifted label back to its tooth. Dodge is per-row (each phase independent),
  // computed left→right over the observed teeth's pixel x (reuses the trace
  // plot's `dodgeX`).
  const observedApexY = y0 - toothH;
  const observed = row.series.teeth
    .map((t, j) => ({ j, px: ctx.x.to(t.q) }))
    .filter((e) => row.series.teeth[e.j]!.observed && Number.isFinite(e.px))
    .sort((a, b) => a.px - b.px);
  const dodged = dodgeX(observed.map((e) => e.px), LABEL_MIN_GAP);
  const labelXByIndex = new Map<number, number>();
  observed.forEach((e, k) => labelXByIndex.set(e.j, dodged[k]!));
  return (
    <g key={i} data-role="comb-row" data-row-kind={row.kind}>
      {row.series.teeth.map((t, j) => {
        const px = ctx.x.to(t.q);
        const hot = isHot(t.q);
        if (t.observed) {
          const apexY = observedApexY;
          const labelX = labelXByIndex.get(j) ?? px;
          const shifted = Math.abs(labelX - px) > 0.5;
          return (
            <g
              key={j}
              data-role="tooth"
              data-q={t.q}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...hov(t.q)}
            >
              <line data-role="tooth-stem" x1={px} y1={y0} x2={px} y2={apexY} stroke={color} strokeWidth={hot ? 2.8 : 2} />
              <circle data-role="tooth-cap" cx={px} cy={apexY} r={hot ? 5 : 2.6} fill={color} />
              {hot ? (
                <circle data-role="tooth-ring" cx={px} cy={apexY} r={8} fill="none" stroke={color} strokeWidth={1.4} opacity={0.6} />
              ) : null}
              {shifted ? (
                <line
                  data-role="tooth-leader"
                  x1={px}
                  y1={apexY - 1}
                  x2={labelX}
                  y2={apexY - 8}
                  stroke="var(--color-hair-strong)"
                  strokeWidth={1}
                />
              ) : null}
              {/* paint-order halo (plate) lifts the ratio off the stems/leaders
                  when the labels are spread over a dense comb. */}
              <text
                data-role="tooth-mlabel" x={labelX} y={apexY - 9} textAnchor="middle"
                className="font-mono" fontSize={9.5} fontWeight={700}
                style={{
                  fill: color,
                  paintOrder: "stroke",
                  stroke: "var(--color-plate)",
                  strokeWidth: 3,
                  strokeLinejoin: "round",
                }}
              >
                {t.label}
              </text>
            </g>
          );
        }
        const ah = toothH * 0.72;
        const apexY = y0 - ah;
        const w = hot ? 2 : 1.3;
        return (
          <g
            key={j}
            data-role="caret"
            data-q={t.q}
            {...(hot ? { "data-hot": "true" } : {})}
            style={{ cursor: onHoverQ ? "pointer" : "default" }}
            {...hov(t.q)}
          >
            <line
              data-role="caret-stem" x1={px} y1={y0} x2={px} y2={apexY}
              stroke={LEFTOVER_COLOR} strokeWidth={w} strokeDasharray="1.5 1.8"
            />
            <path
              data-role="caret-cap" d={`M${px - 3.2} ${apexY} L${px} ${apexY - 4} L${px + 3.2} ${apexY}`}
              fill="none" stroke={LEFTOVER_COLOR} strokeWidth={w}
            />
          </g>
        );
      })}
    </g>
  );
}
