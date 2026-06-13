import { useLayoutEffect, useRef } from "react";
import { type Projection } from "../projection";
import { PeakGlyph } from "../../ui/PeakGlyph";
import { peakGlyph } from "../../ui/peakMark";

/** A one-shot request to move keyboard focus onto a specific peak mark. The
 *  `nonce` lets the consumer re-fire even when the same id repeats (e.g. two
 *  removes in a row that both land focus on the same surviving neighbour).
 *  `id: null` (or an id whose mark no longer exists) means "no peak survived" —
 *  the layer calls `onFocusFallback` so the consumer can park focus elsewhere
 *  (the "+ Peak" toolbar button). */
export interface PeakFocusRequest {
  id: number | null;
  nonce: number;
}

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
  /** One-shot keyboard-focus re-anchor (WCAG 2.4.3). After a destructive edit
   *  unmounts the activated mark, the consumer requests focus on the nearest
   *  surviving peak by id; this layer moves focus to its `<g data-peak-id>` once
   *  per nonce. Only meaningful while editing is armed (focusable marks). */
  focusRequest?: PeakFocusRequest;
  /** Called once per focusRequest nonce when no surviving mark can take focus
   *  (id is null, or its mark is gone). The consumer parks focus elsewhere. */
  onFocusFallback?: () => void;
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
  focusRequest,
  onFocusFallback,
}: PlotPeaksProps): JSX.Element {
  const { x, y } = projection;

  // Consume the focus request once per nonce: after the re-render that dropped
  // the removed mark, move focus onto the surviving peak's <g>. Keyed on the
  // nonce so a foreign/SSE re-render (which leaves the nonce untouched) never
  // re-steals focus — the same wantFocus-ref discipline RepresentativeBox and
  // the Builder RecipeRow use. Layout effect so the move happens before paint.
  const rootRef = useRef<SVGGElement | null>(null);
  const lastNonce = useRef<number | null>(null);
  useLayoutEffect(() => {
    if (!focusRequest) return;
    if (lastNonce.current === focusRequest.nonce) return;
    lastNonce.current = focusRequest.nonce;
    const root = rootRef.current;
    const target =
      root && focusRequest.id != null
        ? root.querySelector<SVGGElement>(
            `g[data-peak-id="${focusRequest.id}"]`,
          )
        : null;
    // A survivor outside the zoom window renders non-focusable (no tabIndex —
    // FO-ZOOMEDIT), so re-anchoring there would silently drop focus to <body>.
    // Fall back to the "+ Peak" button in that case, same as a vanished mark.
    if (target && target.getAttribute("tabindex") === "0") target.focus();
    else onFocusFallback?.();
  }, [focusRequest, onFocusFallback]);

  return (
    <g data-role="plot-peaks" ref={rootRef}>
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
        // FO-ZOOMEDIT: a peak whose q falls outside the visible x-window (after
        // a wheel-zoom) is clipped to invisibility — glyph AND focus ring both
        // gone. Making it a tabbable role=button would seed an invisible tab
        // stop (and pointer can't reach it either, so it would break parity).
        // Detect off-window from the projection's own domain (d3 scales don't
        // clamp, so x.to(q) lands outside the range) and drop its control role.
        const xdom = x.domain;
        const xLo = Math.min(xdom[0], xdom[1]);
        const xHi = Math.max(xdom[0], xdom[1]);
        const inWindow = p.q >= xLo && p.q <= xHi;
        // Hover q-link (onPeakFocus) and keyboard editing (onPeakActivate) are
        // independent gates. Only the latter makes the mark a real control —
        // tabIndex/role/name appear with it and vanish when editing disarms.
        const focusAttrs = onPeakFocus
          ? {
              onFocus: () => onPeakFocus(p.id),
              onBlur: () => onPeakFocus(null),
            }
          : {};
        const activateAttrs = onPeakActivate && inWindow
          ? {
              tabIndex: 0,
              role: "button" as const,
              // FO-RESCORE2 F11: name the peak's PROVENANCE (auto-found vs the
              // user's manual add) and excluded state — it was glyph-only before,
              // invisible to a screen-reader user editing peaks.
              "aria-label": `${p.source === "manual" ? "Manual" : "Auto"} peak at q = ${p.q.toFixed(4)}${p.excluded ? " (excluded)" : ""}`,
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
            // The wrapper carries data-peak-id (mirroring the inner glyph) so
            // the focus-re-anchor layer effect can find the survivor by id.
            data-peak-id={p.id}
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
