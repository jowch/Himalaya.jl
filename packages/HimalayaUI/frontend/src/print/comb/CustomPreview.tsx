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
  /** When set, each observed peak gets a click/Enter hit target that emits its
   *  q — the consumer snaps the lattice so the first reflection lands on it
   *  (click-to-snap). Omit for a static preview. The renderer owns the pointer +
   *  q-geometry; the snap math (lattice-for-first-order) lives in the consumer. */
  onSelectObserved?: (q: number) => void;
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
  onSelectObserved,
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

      {/* observed peak markers (+ optional click-to-snap hit target). The
          visible dot stays 2.5px; the hit circle is a larger transparent
          target so clicking near the dot snaps the first reflection to it. */}
      {observed.map((q, i) => (
        <g key={`obs-${i}`}>
          <circle
            data-observed-q={q}
            cx={x.to(q)}
            cy={baseY}
            r={2.5}
            fill="var(--color-ink-soft)"
          />
          {onSelectObserved && (
            <circle
              data-observed-hit={q}
              cx={x.to(q)}
              cy={baseY}
              r={9}
              fill="transparent"
              role="button"
              tabIndex={0}
              aria-label={`Snap lattice so the first reflection lands on the observed peak at q ${q.toFixed(3)}`}
              style={{ cursor: "pointer" }}
              onClick={() => onSelectObserved(q)}
              onKeyDown={(e) => {
                if (e.key === "Enter" || e.key === " ") {
                  e.preventDefault();
                  onSelectObserved(q);
                }
              }}
            />
          )}
        </g>
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
