import { CombScaffold, ROW_H, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
// Exported so ResidualLegend's band/track note is derived from the SAME constants the
// chart draws — the annotation can't drift from the geometry.
export const RESID_DOMAIN = 0.03; // fixed symmetric y-domain (the "track", ±3%)
export const RESID_BAND = 0.022;  // tolerance band drawn at ±2.2%
// px from baseline to the row's ±RESID_DOMAIN edge. Capped at the room BELOW the
// baseline (baselineY sits ROW_H-14 from the row top → 14px to the bottom edge), so a
// clamped negative residual stays inside the row rather than bleeding into the gap.
const HALF_SPAN = ROW_H / 2 - 10;

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
  // "fit" prefix so the number reads as fit quality, not a bare statistic; the
  // ResidualLegend's band/track note gives the surrounding context.
  const sub = row.series.rSquared !== undefined ? `fit R² ${row.series.rSquared.toFixed(3)}` : undefined;
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
  const bandPx = (RESID_BAND / RESID_DOMAIN) * HALF_SPAN;
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
          const absRes = Math.abs(res);
          // Three tiers:
          //  • in tolerance (within the ±RESID_BAND band)      → filled dot at value
          //  • out of tolerance, still on the track (RESID_BAND..DOMAIN) → HOLLOW dot at value
          //  • off-scale (beyond ±RESID_DOMAIN)                 → chevron clamped at edge, no dot
          // The hollow dot says "out of tolerance" while still reading its real value;
          // the chevron says "ran off the track" (a binary signal — magnitude past the
          // edge is noise on a quality check). Chevron geometry matches the CombChart
          // predicted-absent caret (half-width 3.2, height 4).
          const offScale = absRes > RESID_DOMAIN;
          const outOfTol = !offScale && absRes > RESID_BAND;
          const hot = hoveredQ !== undefined && Math.abs(t.q - hoveredQ) <= TOL;
          const cx = ctx.x.to(t.q);
          if (offScale) {
            const up = res > 0;
            const edgeY = y0 + (up ? -HALF_SPAN : HALF_SPAN);
            const apexY = up ? edgeY - 4 : edgeY + 4; // apex points outward (off the track)
            return (
              <g
                key={j}
                data-role="resid-point"
                data-q={t.q}
                data-overflow="true"
                {...(hot ? { "data-hot": "true" } : {})}
                style={{ cursor: onHoverQ ? "pointer" : "default" }}
                {...(onHoverQ ? { onMouseEnter: () => onHoverQ(t.q), onMouseLeave: () => onHoverQ(undefined) } : {})}
              >
                <path
                  data-role="resid-overflow" data-dir={up ? "up" : "down"}
                  d={`M${cx - 3.2} ${edgeY} L${cx} ${apexY} L${cx + 3.2} ${edgeY}`}
                  fill="none" stroke={color} strokeWidth={hot ? 2 : 1.3}
                  strokeLinejoin="round" strokeLinecap="round"
                />
              </g>
            );
          }
          const cy = yFor(res);
          const r = hot ? 4 : 2.6;
          return (
            <g
              key={j}
              data-role="resid-point"
              data-q={t.q}
              {...(outOfTol ? { "data-outoftol": "true" } : {})}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...(onHoverQ ? { onMouseEnter: () => onHoverQ(t.q), onMouseLeave: () => onHoverQ(undefined) } : {})}
            >
              {/* Hollow when out of tolerance (past the band, still on the track); filled in tolerance. */}
              <circle
                cx={cx} cy={cy} r={r}
                fill={outOfTol ? "none" : color}
                stroke={color} strokeWidth={outOfTol ? 1.6 : hot ? 1.4 : 0}
              />
              {hot ? <circle cx={cx} cy={cy} r={r + 3} fill="none" stroke={color} strokeWidth={1.4} opacity={0.6} /> : null}
            </g>
          );
        })}
    </g>
  );
}
