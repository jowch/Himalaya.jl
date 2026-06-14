import { CombScaffold, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
const LEFTOVER_COLOR = "var(--color-ink-faint)";

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
  const scaffoldRows: ScaffoldRow[] = rows.map(rowToGutter);
  return (
    <div className={className}>
      <CombScaffold
        rows={scaffoldRows}
        xDomain={xDomain}
        ariaLabel="reflection comb"
        {...(maxWidth !== undefined ? { maxWidth } : {})}
      >
        {(ctx) => rows.map((row, i) => renderRow(row, i, ctx, hoveredQ, onHoverQ))}
      </CombScaffold>
    </div>
  );
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
  return (
    <g key={i} data-role="comb-row" data-row-kind={row.kind}>
      {row.series.teeth.map((t, j) => {
        const px = ctx.x.to(t.q);
        const hot = isHot(t.q);
        if (t.observed) {
          const apexY = y0 - toothH;
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
              <text
                data-role="tooth-mlabel" x={px} y={apexY - 7} textAnchor="middle"
                className="font-mono" fontSize={9.5} fontWeight={700} fill={color}
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
