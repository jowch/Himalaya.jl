// peakMark.tsx — <PeakGlyph>, the React/SVG renderer for the peak-glyph
// descriptor (Plan C plot spine). Canonical geometry for BOTH the legend
// (PlotCard) and, via the same path data, the interactive overlay
// (TraceViewer). The §5.1 atoms live in `peakGlyph()` (peakMark.ts); this
// component just paints the descriptor.
//
// Apex orientation is canonical here: the downward triangle points AT the
// peak, the diamond is the manual silhouette, the caret is the hollow
// predicted-absent mark. (Observable Plot's `triangle` symbol points up and
// cannot draw the caret — so on-screen orientation lives in this SVG.)
import type { PeakGlyphDescriptor, PeakMarkOpts } from "./peakMark";
import { peakGlyph } from "./peakMark";

export interface PeakGlyphProps {
  descriptor: PeakGlyphDescriptor;
  /** Centre x in the parent SVG coordinate space. */
  x: number;
  /** Apex / centre y in the parent SVG coordinate space. The glyph is drawn
   *  so its apex (triangle) / centre (diamond/caret) sits at (x, y). */
  y: number;
  /** Stroke used to halo the glyph against the trace (paper by default). The
   *  caller threads the resolved paper colour so the export renderer matches. */
  haloStroke?: string;
  /** Opacity multiplier (ghosting / dim states). */
  opacity?: number;
  /** Optional data attributes for hit-testing / E2E selectors. */
  dataPeakId?: number;
}

/** Build the SVG path geometry for one descriptor centred at (x, y). The
 *  downward triangle's apex is AT y; the diamond/caret are centred on y. */
function glyphPoints(
  shape: PeakGlyphDescriptor["shape"],
  x: number,
  y: number,
  r: number,
): string {
  const hw = r; // half-width
  const h = r * 1.75; // height
  if (shape === "triangle-down") {
    // apex points down (at the peak): two top vertices + the lower apex.
    return `${x - hw},${y - h} ${x + hw},${y - h} ${x},${y}`;
  }
  // diamond: centred on (x, y - h/2) so it floats just above the peak like the
  // triangle's body. Four vertices.
  const cy = y - h / 2;
  return `${x},${cy - hw} ${x + hw},${cy} ${x},${cy + hw} ${x - hw},${cy}`;
}

export function PeakGlyph({
  descriptor,
  x,
  y,
  haloStroke = "var(--color-paper)",
  opacity = 1,
  dataPeakId,
}: PeakGlyphProps): JSX.Element {
  const { shape, fill, stroke, strokeWidth, ring, r } = descriptor;
  const idAttr = dataPeakId !== undefined ? { "data-peak-id": String(dataPeakId) } : {};

  // Hot (q-link) emphasis is carried by the glyph's OWN outline: a darker,
  // heavier accent stroke — NOT a surrounding ring (review note 14).
  return (
    <g data-role="peak-glyph">
      {shape === "caret" ? (
        // Hollow caret: an open chevron pointing down at the predicted q.
        <polyline
          data-shape="caret"
          data-hot={ring ? "true" : undefined}
          points={`${x - r},${y - r * 1.4} ${x},${y} ${x + r},${y - r * 1.4}`}
          fill="none"
          stroke={ring ? "var(--color-accent)" : stroke}
          strokeWidth={ring ? strokeWidth + 1.5 : strokeWidth + 0.5}
          strokeOpacity={opacity}
          strokeLinecap="round"
          strokeLinejoin="round"
          {...idAttr}
        />
      ) : (
        <polygon
          data-shape={shape}
          data-hot={ring ? "true" : undefined}
          points={glyphPoints(shape, x, y, r)}
          fill={fill}
          fillOpacity={fill === "none" ? undefined : opacity}
          stroke={ring ? "var(--color-accent)" : fill === "none" ? stroke : haloStroke}
          strokeWidth={ring ? strokeWidth + 1.5 : strokeWidth}
          strokeOpacity={opacity}
          {...idAttr}
        />
      )}
    </g>
  );
}

/** Convenience: render a glyph straight from atom options (legend call sites). */
export function PeakGlyphFromOpts(
  props: { opts: PeakMarkOpts } & Omit<PeakGlyphProps, "descriptor">,
): JSX.Element {
  const { opts, ...rest } = props;
  return <PeakGlyph descriptor={peakGlyph(opts)} {...rest} />;
}
