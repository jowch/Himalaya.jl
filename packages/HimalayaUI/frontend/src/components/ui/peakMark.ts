// peakMark.ts — the single peak-glyph builder (Plan C plot spine).
//
// One source of truth for the §5.1 peak-encoding atoms. It produces BOTH:
//   - `peakGlyph(opts)` → a pure geometry descriptor consumed by the SVG
//     overlay (`<PeakGlyph>`) AND the legend. The descriptor is canonical for
//     on-screen apex orientation (downward triangle, diamond, caret).
//   - `peakMark(rows)` → an Observable Plot `Markish` (Plot.dot) for the
//     export / member-thumbnail contexts where exact apex orientation is
//     non-critical. Caret/predicted-absent glyphs MUST go through the SVG
//     `peakGlyph` path — `Plot.dot` cannot draw the caret.
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
// COLOUR CONTRACT: callers pass a RESOLVED colour string (phaseColor(),
// LIGHT_PALETTE.*, or a CSS var token) — peakMark NEVER reads a CSS var here,
// so the export/canvas renderer (no CSS vars) and the on-screen overlay share
// one builder. Magenta `--color-peak-manual` is RETIRED: manual provenance is
// the diamond silhouette, not a hue.
import * as Plot from "@observablehq/plot";

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
  /** terracotta hot ring (q-link). */
  ring: boolean;
  r: number;
  offsetPx: number;
  /** false when optimistic — the overlay must not route clicks here. */
  interactive: boolean;
}

const HOT_GROWTH = 1.5;
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
  const baseR = opts.r ?? DEFAULT_R;
  const r = opts.hot ? baseR * HOT_GROWTH : baseR;
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

/** A peak row for the Observable Plot Markish. Each row carries q, y and a
 *  RESOLVED colour; state flags are optional. */
export type PeakMarkRow = { q: number; y: number; color: string } & PeakMarkState;

/** Observable Plot Markish for a set of peak rows (export + member layers).
 *
 *  NOTE: Plot has no native downward-triangle symbol matching the overlay
 *  geometry, and Plot.dot cannot draw the predicted-absent CARET. So this
 *  markish is used ONLY where exact apex orientation is non-critical
 *  (export / member thumbnails). For caret atoms and on-screen orientation,
 *  the SVG `peakGlyph`/`<PeakGlyph>` path is canonical. */
export function peakMark(
  rows: ReadonlyArray<PeakMarkRow>,
  channels: { x?: string; y?: string; r?: number } = {},
): Plot.Markish {
  return Plot.dot(rows as unknown as Plot.Data, {
    x: channels.x ?? "q",
    y: channels.y ?? "y",
    r: channels.r ?? DEFAULT_R,
    symbol: (d: unknown) =>
      (d as PeakMarkState).source === "manual" ? "diamond" : "triangle",
    fill: (d: unknown) => {
      const s = d as PeakMarkRow;
      return s.predictedAbsent || s.optimistic || s.excluded ? "none" : s.color;
    },
    stroke: (d: unknown) => (d as PeakMarkRow).color,
    strokeWidth: 0.75,
    fillOpacity: (d: unknown) => ((d as PeakMarkState).excluded ? 0.35 : 1),
  });
}
