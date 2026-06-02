import { phaseColor } from "../../phases";
import { makeAxis, positiveExtent } from "../plot/projection";
import type { CombSeries } from "./combModel";

export interface CustomPreviewProps {
  /** Candidate phase's predicted comb at the chosen lattice. */
  series: CombSeries;
  /** Observed peak q-values (Å⁻¹) to overlay. */
  observed: number[];
  /** SVG width (px). */
  width?: number;
  /** SVG height (px). */
  height?: number;
  /** Placement only. */
  className?: string;
}

const mL = 10;
const mR = 10;
const mT = 14;
const mB = 26;

/**
 * A static single-row comb showing one candidate phase's predicted Bragg teeth
 * (at a chosen lattice) overlaid against the observed peak q-values. The live
 * preview inside the custom-index modal — hand-rolled (no gutter/scroll/hover)
 * so it stays light; reuses `makeAxis` for the log-q scale and `phaseColor()`
 * for the phase hue.
 */
export function CustomPreview({
  series,
  observed,
  width = 520,
  height = 150,
  className,
}: CustomPreviewProps): JSX.Element {
  const pw = width - mL - mR;
  const baseY = height - mB;
  const topY = mT;

  // Padded log-q domain over every predicted tooth + observed peak (~8% each side).
  const qs = [...series.teeth.map((t) => t.q), ...observed];
  const [lo, hi] = positiveExtent(qs, [0.01, 0.2]);
  const x = makeAxis([lo * 0.92, hi * 1.08], [mL, mL + pw], "log");

  const color = phaseColor(series.phase);

  return (
    <svg
      viewBox={`0 0 ${width} ${height}`}
      role="img"
      aria-label="Custom index reflections against observed peaks"
      data-testid="custom-preview"
      data-phase={series.phase}
      className={className}
    >
      {/* baseline */}
      <line
        x1={mL}
        y1={baseY}
        x2={mL + pw}
        y2={baseY}
        stroke="var(--color-hair-strong)"
        strokeWidth={1}
      />

      {/* observed peak markers */}
      {observed.map((q, i) => (
        <circle
          key={`obs-${i}`}
          data-observed-q={q}
          cx={x.to(q)}
          cy={baseY}
          r={2.5}
          fill="var(--color-ink-soft)"
        />
      ))}

      {/* predicted teeth */}
      {series.teeth.map((t, i) => {
        const px = x.to(t.q);
        if (t.observed) {
          return (
            <g key={`tooth-${i}`}>
              <line
                data-tooth-q={t.q}
                data-tooth-observed={t.observed}
                data-tooth-label={t.label}
                x1={px}
                y1={baseY}
                x2={px}
                y2={topY}
                style={{ stroke: color }}
                strokeWidth={1.8}
                strokeLinecap="round"
              />
              <text
                x={px}
                y={topY - 3}
                textAnchor="middle"
                fontFamily="var(--font-mono)"
                fontSize={8}
                fill="var(--color-ink-faint)"
              >
                {t.label}
              </text>
            </g>
          );
        }
        return (
          <line
            key={`caret-${i}`}
            data-tooth-q={t.q}
            data-tooth-observed={t.observed}
            data-tooth-label={t.label}
            x1={px}
            y1={baseY}
            x2={px}
            y2={topY}
            stroke="var(--color-ink-faint)"
            strokeDasharray="1.5 1.8"
            strokeWidth={1.5}
          />
        );
      })}
    </svg>
  );
}
