import { type Projection } from "../projection";
import { PeakGlyph } from "../../ui/PeakGlyph";
import { peakGlyph } from "../../ui/peakMark";

export interface PlotPeak {
  id: number;
  q: number;
  intensity?: number | null;
  source: "auto" | "manual";
  excluded?: boolean;
  predictedAbsent?: boolean;
  hot?: boolean;
  /** Per-peak resolved colour (e.g. assigned-index colour). Falls back to the layer `color` prop. */
  color?: string;
  /** Pre-resolved label string (consumed by label layers; unused in this mark). */
  label?: string;
}

export interface PlotPeaksProps {
  peaks: PlotPeak[];
  projection: Projection;
  /** Resolved colour for the glyphs (phase colour, usually). */
  color: string;
  /** y data-value for peaks lacking an intensity (anchors them to baseline). */
  baselineI?: number;
  /** Paper colour threaded to PeakGlyph's halo (export-parity). */
  paperColor?: string;
  /** Keyboard focus/blur handler for accessibility. Called with peak id on focus, null on blur. */
  onPeakFocus?: (id: number | null) => void;
  /** Armed keyboard editing. Its PRESENCE gates focusability: each mark becomes
   *  a tabbable role=button named by its q. Enter/Space activate with
   *  altKey=false (remove); Alt+Enter activates with altKey=true (exclude).
   *  Omit it while editing is disarmed so read-only marks stay roleless — an
   *  inert focusable control would lie. */
  onPeakActivate?: (id: number, altKey: boolean) => void;
  /** When non-empty, peaks NOT in this set (and not hot) fade to neutral gray. */
  highlightPeakIds?: ReadonlySet<number>;
}

export function PlotPeaks({
  peaks,
  projection,
  color,
  baselineI,
  paperColor,
  onPeakFocus,
  onPeakActivate,
  highlightPeakIds,
}: PlotPeaksProps): JSX.Element {
  const { x, y } = projection;
  return (
    <g data-role="plot-peaks">
      {peaks.map((p) => {
        if (p.id < 0) return null;
        const px = x.to(p.q);
        if (!Number.isFinite(px)) return null;
        const iVal = p.intensity ?? baselineI;
        const py =
          iVal != null && Number.isFinite(iVal) ? y.to(iVal) : y.range[0];
        const hl = highlightPeakIds;
        const dimmed = !!hl && hl.size > 0 && !hl.has(p.id) && !p.hot; // hot wins over dim
        const c = dimmed ? "var(--color-ink-faint)" : (p.color ?? color);
        // The resolved colour lands on the wrapping <g>'s CSS `color`; the
        // glyph + q-line paint via currentColor. SVG fill/stroke attributes
        // are NOT in the global transition-property list, but `color` IS
        // (styles.css @layer base: 120ms ease-out) — so routing the paint
        // through `color` makes the losing-peak dim (and its release) ease
        // instead of snapping. prefers-reduced-motion is handled by the
        // global near-zero rule; the data-dimmed flip itself stays instant.
        const descriptor = peakGlyph({
          source: p.source,
          color: "currentColor",
          ...(p.predictedAbsent ? { predictedAbsent: true } : {}),
          ...(p.excluded ? { excluded: true } : {}),
          ...(p.hot ? { hot: true } : {}),
        });
        // Hover q-link (onPeakFocus) and keyboard editing (onPeakActivate) are
        // independent gates. Only the latter makes the mark a real control —
        // tabIndex/role/name appear with it and vanish when editing disarms.
        const focusAttrs = onPeakFocus
          ? {
              onFocus: () => onPeakFocus(p.id),
              onBlur: () => onPeakFocus(null),
            }
          : {};
        const activateAttrs = onPeakActivate
          ? {
              tabIndex: 0,
              role: "button" as const,
              "aria-label": `Peak at q = ${p.q.toFixed(4)}`,
              "aria-keyshortcuts": "Enter Space Alt+Enter",
              onKeyDown: (e: React.KeyboardEvent<SVGGElement>) => {
                if (e.key !== "Enter" && e.key !== " ") return;
                // Holding the key must not machine-gun mutations.
                if (e.repeat) return;
                // preventDefault so Space doesn't scroll the page.
                e.preventDefault();
                // Alt modifies Enter only — the advertised shortcuts are
                // "Enter Space Alt+Enter"; Alt+Space stays a plain remove.
                onPeakActivate(p.id, e.altKey && e.key === "Enter");
              },
            }
          : {};
        return (
          <g
            key={p.id}
            {...focusAttrs}
            {...activateAttrs}
            {...(dimmed ? { "data-dimmed": "true" } : {})}
            style={{ color: c }}
          >
            {p.hot ? (
              // Drop a guide DOWN from the marker to the axis baseline (where
              // the q-readout chip sits) — it no longer runs through the mark.
              <line
                data-role="peak-qline"
                x1={px}
                y1={py}
                x2={px}
                y2={y.range[0]}
                stroke="currentColor"
                strokeWidth={1.5}
                opacity={0.6}
              />
            ) : null}
            <PeakGlyph
              descriptor={descriptor}
              x={px}
              y={py}
              dataPeakId={p.id}
              {...(paperColor ? { haloStroke: paperColor } : {})}
            />
          </g>
        );
      })}
    </g>
  );
}
