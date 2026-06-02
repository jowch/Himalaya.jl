import { CombScaffold, ROW_H, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
const RESID_DOMAIN = 0.03; // fixed symmetric y-domain (±3%)
const BAND = 0.022;        // tolerance band drawn at ±2.2%
const HALF_SPAN = ROW_H / 2 - 9; // px from baseline to the row's ±RESID_DOMAIN edge

interface Props {
  assigned: CombSeries[];
  hovered?: CombSeries;
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  maxWidth?: number;
  className?: string;
}

export function ResidualChart({ assigned, hovered, hoveredQ, onHoverQ, maxWidth, className }: Props): JSX.Element {
  // No leftover row in the residual view (unindexed peaks have no predicted q).
  const rows = assembleRows(assigned, hovered, []);
  const xDomain = combQDomain(rows);
  const scaffoldRows: ScaffoldRow[] = rows.map(rowToGutter);
  return (
    <div className={className}>
      <CombScaffold
        rows={scaffoldRows}
        xDomain={xDomain}
        ariaLabel="indexing residuals"
        {...(maxWidth !== undefined ? { maxWidth } : {})}
      >
        {(ctx) => rows.map((row, i) => renderResidRow(row, i, ctx, hoveredQ, onHoverQ))}
      </CombScaffold>
    </div>
  );
}

function rowToGutter(row: CombRow): ScaffoldRow {
  if (row.kind === "leftover") return { gutterTitle: "leftover" }; // unreachable (no leftover passed)
  const sub = row.series.rSquared !== undefined ? `R² ${row.series.rSquared.toFixed(3)}` : undefined;
  const base: ScaffoldRow = { gutterTitle: row.series.phase, ...(sub ? { gutterSub: sub } : {}) };
  return row.kind === "preview" ? { ...base, preview: true } : base;
}

function renderResidRow(
  row: CombRow,
  i: number,
  ctx: ScaffoldCtx,
  hoveredQ: number | undefined,
  onHoverQ?: (q?: number) => void,
): JSX.Element {
  if (row.kind === "leftover") return <g key={i} data-role="resid-row" data-row-kind="leftover" />; // unreachable
  const { color } = row.series;
  const y0 = ctx.baselineY(i);
  const bandPx = (BAND / RESID_DOMAIN) * HALF_SPAN;
  // Positive residual (over-predicted) sits ABOVE the baseline (smaller y).
  const yFor = (res: number): number => {
    const clamped = Math.max(-RESID_DOMAIN, Math.min(RESID_DOMAIN, res));
    return y0 - (clamped / RESID_DOMAIN) * HALF_SPAN;
  };
  return (
    <g key={i} data-role="resid-row" data-row-kind={row.kind}>
      <rect
        data-role="tolerance-band"
        x={ctx.padX} y={y0 - bandPx} width={ctx.plotW} height={bandPx * 2}
        fill="var(--color-ink)" opacity={0.05}
      />
      {row.series.teeth
        .filter((t) => t.observed && t.residual !== undefined)
        .map((t, j) => {
          const res = t.residual as number;
          const overflow = Math.abs(res) > RESID_DOMAIN;
          const hot = hoveredQ !== undefined && Math.abs(t.q - hoveredQ) <= TOL;
          const cx = ctx.x.to(t.q);
          const cy = yFor(res);
          const r = hot ? 4 : 2.6;
          return (
            <g
              key={j}
              data-role="resid-point"
              data-q={t.q}
              {...(overflow ? { "data-overflow": "true" } : {})}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...(onHoverQ ? { onMouseEnter: () => onHoverQ(t.q), onMouseLeave: () => onHoverQ(undefined) } : {})}
            >
              <circle
                cx={cx} cy={cy} r={r}
                fill={overflow ? "none" : color}
                stroke={color}
                strokeWidth={overflow ? 1.6 : hot ? 1.4 : 0}
              />
              {hot ? <circle cx={cx} cy={cy} r={r + 3} fill="none" stroke={color} strokeWidth={1.4} opacity={0.6} /> : null}
            </g>
          );
        })}
    </g>
  );
}
