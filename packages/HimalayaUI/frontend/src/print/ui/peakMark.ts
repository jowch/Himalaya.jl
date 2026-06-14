// peakMark.ts — the single peak-glyph builder (Plan C plot spine).
//
// One source of truth for the §5.1 peak-encoding atoms: `peakGlyph(opts)` → a
// pure geometry descriptor consumed by the SVG overlay (`<PeakGlyph>`) AND the
// legend. The descriptor is canonical for on-screen apex orientation (downward
// triangle, diamond, caret).
//
// §5.1 encoding atoms (see docs plan):
//   colour = phase, ALWAYS. Provenance/state via silhouette + fill, never hue.
//   auto              → filled DOWNWARD triangle (points at the peak)
//   manual            → DIAMOND (neutral-gray when unindexed; phase once indexed)
//   excluded          → ghosted (hollow), not struck
//   predicted-absent  → hollow CARET
//   optimistic (id<0) → outline-only (no fill), non-interactive
//   hot (q-link)      → grow + ring, NOT a recolour
//
// COLOUR CONTRACT: callers pass a RESOLVED colour string (phaseColor() or a CSS
// var token) — peakGlyph NEVER reads a CSS var here, so the canvas renderer (no
// CSS vars) and the on-screen overlay share one builder. Magenta
// `--color-peak-manual` is RETIRED: manual provenance is the diamond
// silhouette, not a hue.

export type PeakSource = "auto" | "manual";

export interface PeakMarkState {
  source: PeakSource;
  /** ghosted (hollow, faint) — user-excluded auto peak. */
  excluded?: boolean;
  /** id<0 → outline-only, non-interactive (pending server confirmation). */
  optimistic?: boolean;
  /** hollow caret — a predicted q with no observed peak. */
  predictedAbsent?: boolean;
  /** q-link: grow + ring (no recolour). */
  hot?: boolean;
}

export interface PeakMarkOpts extends PeakMarkState {
  /** Resolved colour: phase colour when indexed, neutral-gray when unindexed
   *  manual. NEVER a CSS-var READ here — callers pass the resolved string. */
  color: string;
  /** Pixel offset of the glyph apex above the baseline/peak. 7px standalone,
   *  5px member — parameterized so the existing feel doesn't re-drift. */
  offsetPx?: number;
  /** Base glyph radius in px (grows when hot). */
  r?: number;
}

export interface PeakGlyphDescriptor {
  shape: "triangle-down" | "diamond" | "caret";
  fill: string | "none";
  stroke: string;
  strokeWidth: number;
  /** q-link emphasis flag — the q-line/readout carry it; the mark is unchanged
   *  (kept as a `data-hot` DOM hook). */
  ring: boolean;
  r: number;
  offsetPx: number;
  /** false when optimistic — the overlay must not route clicks here. */
  interactive: boolean;
}

const DEFAULT_R = 4;
const DEFAULT_OFFSET_PX = 7;

/** Resolve the §5.1 atoms into a geometry descriptor — the single source of
 *  truth for the SVG overlay AND the legend. */
export function peakGlyph(opts: PeakMarkOpts): PeakGlyphDescriptor {
  const shape: PeakGlyphDescriptor["shape"] = opts.predictedAbsent
    ? "caret"
    : opts.source === "manual"
      ? "diamond"
      : "triangle-down";
  // Hot does NOT change the mark at all: the q-link emphasis is the q-line +
  // q-readout (PlotPeaks/TracePlot), so the glyph keeps its resting size+colour.
  const r = opts.r ?? DEFAULT_R;
  // Filled only when it's a real, present, non-ghosted, non-optimistic peak.
  const filled = !opts.predictedAbsent && !opts.optimistic && !opts.excluded;
  return {
    shape,
    fill: filled ? opts.color : "none",
    stroke: opts.color,
    // Outline-only atoms (optimistic) draw a heavier stroke so the empty glyph
    // still reads against the trace.
    strokeWidth: opts.optimistic ? 1.25 : 0.75,
    ring: opts.hot === true,
    r,
    offsetPx: opts.offsetPx ?? DEFAULT_OFFSET_PX,
    interactive: opts.optimistic !== true,
  };
}

// Re-export the SVG renderer so consumers can import the whole peak-glyph
// surface (builder + descriptor + component) from one module specifier. The
// component lives in peakMark.tsx (JSX); this binding keeps `./peakMark` the
// single import path. (ESM tolerates the .ts↔.tsx cycle — .tsx only consumes
// the type + peakGlyph above, both fully initialized before the component runs.)
export { PeakGlyph, PeakGlyphFromOpts } from "./PeakGlyph";
