import { type Projection } from "../projection";
import { type PlotPeak } from "./PlotPeaks";

export interface PlotLabelsProps {
  peaks: PlotPeak[];
  projection: Projection;
  /** Fallback colour when a peak has no per-peak colour. */
  color: string;
  /** Pixel offset above the peak y-position where the label renders. */
  offsetPx?: number;
  /** Minimum horizontal pixel gap enforced between adjacent labels. */
  labelWidthPx?: number;
  /** y data-value for peaks lacking an intensity (anchors them to baseline). */
  baselineI?: number;
  /** When non-empty, labels for peaks NOT in this set fade to neutral gray. */
  highlightPeakIds?: ReadonlySet<number>;
}

/**
 * Pure two-pass horizontal dodge for a pre-sorted ascending array of pixel
 * x-positions.
 *
 * Pass 1 (left→right sweep): each position must be at least `labelWidthPx`
 * to the right of the previous. This guarantees no two labels overlap but
 * shifts the whole cluster rightward.
 *
 * Pass 2 (recenter): shift all dodged positions by (meanNatural − meanDodged)
 * so the cluster stays symmetric about the original mean. Gaps are preserved
 * because a uniform translation keeps all inter-element distances the same.
 */
export function dodgeX(xs: number[], labelWidthPx: number): number[] {
  if (xs.length === 0) return [];
  if (xs.length === 1) return [xs[0]!];

  // Pass 1: left-to-right sweep.
  const out: number[] = new Array(xs.length);
  out[0] = xs[0]!;
  for (let i = 1; i < xs.length; i++) {
    out[i] = Math.max(xs[i]!, out[i - 1]! + labelWidthPx);
  }

  // Pass 2: recenter so mean(dodged) === mean(natural).
  let naturalSum = 0;
  let dodgedSum = 0;
  for (let i = 0; i < xs.length; i++) {
    naturalSum += xs[i]!;
    dodgedSum += out[i]!;
  }
  const shift = (naturalSum - dodgedSum) / xs.length;
  if (shift !== 0) {
    for (let i = 0; i < xs.length; i++) {
      out[i] = out[i]! + shift;
    }
  }

  return out;
}

/** One renderable label entry (pixel-space, pre-sorted). */
interface LabelEntry {
  id: number;
  naturalPx: number;
  py: number;
  label: string;
  color: string;
  dimmed: boolean;
}

export function PlotLabels({
  peaks,
  projection,
  color,
  offsetPx = 12,
  labelWidthPx = 30,
  baselineI,
  highlightPeakIds,
}: PlotLabelsProps): JSX.Element {
  const { x, y } = projection;

  // Filter to renderable peaks, then map to pixel entries.
  const entries: LabelEntry[] = [];
  for (const p of peaks) {
    if (p.id < 0) continue;
    if (p.excluded) continue;
    if (p.label == null || p.label === "") continue;

    const px = x.to(p.q);
    if (!Number.isFinite(px)) continue;

    const iVal = p.intensity ?? baselineI;
    const py =
      iVal != null && Number.isFinite(iVal) ? y.to(iVal) : y.range[0];

    const dimmed =
      !!highlightPeakIds &&
      highlightPeakIds.size > 0 &&
      !highlightPeakIds.has(p.id);
    entries.push({
      id: p.id,
      naturalPx: px,
      py,
      label: p.label,
      color: dimmed ? "var(--color-ink-faint)" : (p.color ?? color),
      dimmed,
    });
  }

  if (entries.length === 0) {
    return <g data-role="plot-labels" />;
  }

  // Sort ascending by natural pixel x (required by dodgeX).
  entries.sort((a, b) => a.naturalPx - b.naturalPx);

  const naturalXs = entries.map((e) => e.naturalPx);
  const dodgedXs = dodgeX(naturalXs, labelWidthPx);

  return (
    <g data-role="plot-labels">
      {entries.map((e, i) => {
        const dodgedX = dodgedXs[i]!;
        const naturalX = e.naturalPx;
        const labelY = e.py - offsetPx;
        const shifted = Math.abs(dodgedX - naturalX) > 0.5;

        return (
          <g key={e.id}>
            {shifted ? (
              <line
                x1={naturalX}
                y1={e.py - 6}
                x2={dodgedX}
                y2={e.py - offsetPx + 2}
                stroke="var(--color-hair)"
                strokeWidth={1}
              />
            ) : null}
            <text
              data-role="peak-label"
              x={dodgedX}
              y={labelY}
              textAnchor="middle"
              {...(e.dimmed ? { "data-dimmed": "true" } : {})}
              style={{
                fontFamily: "var(--font-mono)",
                fontSize: 11,
                fontWeight: 700,
                fill: e.color,
              }}
            >
              {e.label}
            </text>
          </g>
        );
      })}
    </g>
  );
}
