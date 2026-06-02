import type { WaterfallRow } from "../waterfall/waterfallModel";
import { makeAxis, axisTicks, positiveExtent } from "../plot/projection";

export interface CleanFigureProps {
  rows: WaterfallRow[];
  title: string;
  footer: string;
  /** SVG viewBox width. */
  width?: number;
  /** SVG viewBox height. */
  height?: number;
  /** Placement-only className on the outer card. */
  className?: string;
}

const ARIAL = "Arial, Helvetica, sans-serif";

// Pure phase colours in plain hex — the un-branded export idiom (NOT Print OKLCH).
const HEX: Record<string, string> = {
  Pn3m: "#c8841f",
  Lamellar: "#4a4ba8",
  Im3m: "#4a4ba8",
};

/**
 * The plain GraphPad-style EXPORT renderer. Deliberately sheds "The Print"
 * branding: Arial, literal hex colours, black axes — a self-contained
 * scientific figure for a slide or a colleague.
 */
export function CleanFigure({
  rows,
  title,
  footer,
  width = 600,
  height = 340,
  className,
}: CleanFigureProps) {
  const m = { l: 52, r: 64, t: 12, b: 40 };
  const pw = width - m.l - m.r;
  const ph = height - m.t - m.b;
  const baseY = m.t + ph;

  // Log-q x-axis over all rows' q, padded ~8% each side.
  const [qLo, qHi] = positiveExtent(rows.flatMap((r) => r.trace.q));
  const padded: [number, number] = [qLo * 0.92, qHi * 1.08];
  const x = makeAxis(padded, [m.l, m.l + pw], "log");
  const inDomain = (q: number) => q >= padded[0] && q <= padded[1];
  const ticks = axisTicks(x);

  // Global max intensity across ALL rows, for shared vertical normalization.
  let globalMaxI = 0;
  for (const r of rows) {
    for (const v of r.trace.I) {
      if (Number.isFinite(v) && v > globalMaxI) globalMaxI = v;
    }
  }
  const safeMaxI = globalMaxI > 0 ? globalMaxI : 1;
  const slot = rows.length > 0 ? ph / rows.length : ph;
  const amp = slot * 0.9;

  const pathFor = (row: WaterfallRow, rowBase: number): string => {
    const { q, I } = row.trace;
    let d = "";
    let started = false;
    for (let k = 0; k < q.length; k++) {
      const qk = q[k]!;
      if (!inDomain(qk)) continue;
      const normI = (I[k] ?? 0) / safeMaxI;
      const px = x.to(qk);
      const py = rowBase - normI * amp;
      d += `${started ? "L" : "M"}${px.toFixed(1)} ${py.toFixed(1)} `;
      started = true;
    }
    return d.trim();
  };

  return (
    <div
      data-testid="clean-figure"
      className={className}
      style={{
        background: "#ffffff",
        border: "1px solid #e3e3e3",
        borderRadius: 3,
        padding: "16px 14px 10px",
      }}
    >
      <div
        style={{
          font: `700 14px ${ARIAL}`,
          color: "#111",
          textAlign: "center",
        }}
      >
        {title}
      </div>

      <svg
        viewBox={`0 0 ${width} ${height}`}
        role="img"
        aria-label={title}
        style={{ display: "block", width: "100%" }}
      >
        {/* Black L-shaped axes — GraphPad idiom. */}
        <line x1={m.l} y1={baseY} x2={m.l + pw} y2={baseY} stroke="#111" strokeWidth={1.6} />
        <line x1={m.l} y1={m.t} x2={m.l} y2={baseY} stroke="#111" strokeWidth={1.6} />

        {/* X ticks + labels (labels only for non-minor ticks). */}
        {ticks.map((t, i) => {
          if (!inDomain(t.value)) return null;
          const px = x.to(t.value);
          return (
            <g key={`xt-${i}`}>
              <line x1={px} y1={baseY} x2={px} y2={baseY + 5} stroke="#111" strokeWidth={1} />
              {t.kind !== "minor" && (
                <text
                  x={px}
                  y={baseY + 16}
                  textAnchor="middle"
                  fontFamily={ARIAL}
                  fontSize={9}
                  fill="#111"
                >
                  {t.value.toFixed(2)}
                </text>
              )}
            </g>
          );
        })}

        {/* X axis title. */}
        <text
          x={m.l + pw / 2}
          y={height - 4}
          textAnchor="middle"
          fontFamily={ARIAL}
          fontSize={11}
          fontWeight={700}
          fill="#111"
        >
          {"q (Å⁻¹)"}
        </text>

        {/* Y axis title, rotated. */}
        <text
          x={13}
          y={m.t + ph / 2}
          textAnchor="middle"
          transform={`rotate(-90 13 ${m.t + ph / 2})`}
          fontFamily={ARIAL}
          fontSize={11}
          fontWeight={700}
          fill="#111"
        >
          Intensity (a.u.)
        </text>

        {/* Traces, stacked bottom→top. */}
        {rows.map((row, i) => {
          const rowBase = baseY - i * slot;
          const d = pathFor(row, rowBase);
          return (
            <g key={row.key}>
              <path
                d={d}
                fill="none"
                stroke={(row.phase != null && HEX[row.phase]) || "#777"}
                strokeWidth={2}
                data-phase={row.phase ?? "none"}
              />
              <text
                x={m.l + pw + 6}
                y={rowBase}
                fontFamily={ARIAL}
                fontSize={9}
                fill="#111"
              >
                {row.label}
              </text>
            </g>
          );
        })}
      </svg>

      <div
        style={{
          font: `9px ${ARIAL}`,
          color: "#666",
          textAlign: "center",
        }}
      >
        {footer}
      </div>
    </div>
  );
}
