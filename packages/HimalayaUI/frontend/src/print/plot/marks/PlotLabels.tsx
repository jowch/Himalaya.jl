import { type Projection } from "../projection";
import { type PlotPeak } from "./PlotPeaks";
import { measureTextWidth } from "../../../lib/plot/measureText";

// The label type face (mirrors the <text> style below) — fed to the text
// measurer so the dodge knows each label's real width.
const LABEL_FONT = { px: 11, weight: 700, family: "monospace" } as const;
/** Minimum clear gap between adjacent label edges before the dodge separates them. */
const LABEL_GAP_PX = 5;
/** The marker's vertical centre sits this far above the peak's y (MARKER_LIFT 3
 *  + half the glyph height ~3.5). The connector attaches HERE — the marker's
 *  middle — not at its downward apex ("bottom point"). */
const MARKER_CENTER_DY = 6.5;

export interface PlotLabelsProps {
  peaks: PlotPeak[];
  projection: Projection;
  /** Fallback colour when a peak has no per-peak colour. */
  color: string;
  /** Pixel offset above the peak y-position where the label renders. */
  offsetPx?: number;
  /** y data-value for peaks lacking an intensity (anchors them to baseline). */
  baselineI?: number;
  /** When non-empty, labels for peaks NOT in this set fade to neutral gray. */
  highlightPeakIds?: ReadonlySet<number>;
  /** Z-order split so the dodge connector can render BENEATH the peak markers
   *  (TracePlot draws `"connectors"` before PlotPeaks and `"text"` after).
   *  Default `"all"` renders both in one pass. */
  part?: "all" | "connectors" | "text";
}

/**
 * Pure horizontal dodge for a pre-sorted ascending array of pixel x-centers,
 * width-aware so it only moves labels that ACTUALLY overlap.
 *
 * `halfWidths[i]` is half of label i's rendered width; `gapPx` is the minimum
 * clear gap between adjacent label EDGES. A single forward sweep pushes each
 * label just past the previous label's right edge when — and only when — they
 * would collide. A well-spaced label keeps its natural x exactly, so it draws
 * no spurious connector.
 *
 * (The previous version used a fixed worst-case width AND a global recenter,
 * which translated EVERY label off its natural x — so even cleanly-separated
 * labels registered as "dodged" and grew a connector. That was the aggressive
 * dodging this replaces.)
 */
export function dodgeX(
  xs: number[],
  halfWidths: number[],
  gapPx: number,
): number[] {
  const n = xs.length;
  if (n === 0) return [];
  const out = xs.slice();
  for (let i = 1; i < n; i++) {
    const minCenter = out[i - 1]! + halfWidths[i - 1]! + halfWidths[i]! + gapPx;
    if (out[i]! < minCenter) out[i] = minCenter;
  }
  return out;
}

/** One renderable label entry (pixel-space, pre-sorted). */
interface LabelEntry {
  id: number;
  naturalPx: number;
  py: number;
  label: string;
  /** Half the label's measured rendered width (for width-aware dodging). */
  halfWidth: number;
  color: string;
  dimmed: boolean;
}

export function PlotLabels({
  peaks,
  projection,
  color,
  offsetPx = 12,
  baselineI,
  highlightPeakIds,
  part = "all",
}: PlotLabelsProps): JSX.Element {
  const { x, y } = projection;

  // Filter to renderable peaks, then map to pixel entries (with measured width
  // so the dodge only separates labels that genuinely overlap).
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
      halfWidth: measureTextWidth(p.label, LABEL_FONT) / 2,
      color: dimmed ? "var(--color-ink-faint)" : (p.color ?? color),
      dimmed,
    });
  }

  if (entries.length === 0) {
    return <g data-role="plot-labels" />;
  }

  // Sort ascending by natural pixel x (required by dodgeX).
  entries.sort((a, b) => a.naturalPx - b.naturalPx);

  const dodgedXs = dodgeX(
    entries.map((e) => e.naturalPx),
    entries.map((e) => e.halfWidth),
    LABEL_GAP_PX,
  );

  return (
    <g data-role="plot-labels">
      {entries.map((e, i) => {
        const dodgedX = dodgedXs[i]!;
        const naturalX = e.naturalPx;
        const labelY = e.py - offsetPx;
        const shifted = Math.abs(dodgedX - naturalX) > 0.5;

        return (
          <g key={e.id}>
            {/* Connector only when the label actually moved. It attaches at the
                marker's CENTRE (naturalX, py − MARKER_CENTER_DY) and renders in
                the `"connectors"` pass, which TracePlot draws BEFORE the markers
                — so the marker sits on top and the leader tucks underneath it. */}
            {part !== "text" && shifted ? (
              <line
                x1={naturalX}
                y1={e.py - MARKER_CENTER_DY}
                x2={dodgedX}
                y2={e.py - offsetPx + 2}
                stroke="var(--color-hair)"
                strokeWidth={1}
              />
            ) : null}
            {part === "connectors" ? null : (
            <text
              data-role="peak-label"
              x={dodgedX}
              y={labelY}
              textAnchor="middle"
              {...(e.dimmed ? { "data-dimmed": "true" } : {})}
              // The resolved colour lands on CSS `color`, painted via
              // currentColor: `color` is in the global 120ms ease-out
              // transition list (styles.css @layer base) while `fill` is not,
              // so the dim (and its release) eases instead of snapping —
              // matching the PlotPeaks glyph treatment. Reduced motion is
              // handled by the global near-zero rule.
              // BU-PEAKORD-HALO: a plate-colour halo (stroke painted UNDER the
              // fill via paint-order) lifts the ordinal off the gridlines on
              // low-amplitude traces, where it would otherwise merge with them.
              style={{
                fontFamily: "var(--font-mono)",
                fontSize: 11,
                fontWeight: 700,
                color: e.color,
                fill: "currentColor",
                paintOrder: "stroke",
                stroke: "var(--color-plate)",
                strokeWidth: 3,
                strokeLinejoin: "round",
              }}
            >
              {e.label}
            </text>
            )}
          </g>
        );
      })}
    </g>
  );
}
