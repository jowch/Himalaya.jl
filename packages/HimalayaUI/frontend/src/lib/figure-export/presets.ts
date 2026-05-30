// presets.ts — pinned export dimensions, palette, fonts, strokes.
// Numeric values live here so adapters stay terse and figure-size changes
// are a one-line edit.

export const TRACE_DIMS = { width: 800, height: 600 } as const;
export const COMPARE_DIMS = { width: 1000, height: 400 } as const;

/** Margin defaults applied by the renderer (overridable via spec.plot). */
export const EXPORT_MARGIN = {
  top: 60,    // room for primary + secondary title text
  right: 24,
  bottom: 60, // room for legend rows below the plot
  left: 60,
} as const;

/** Plot title font metrics. */
export const TITLE_FONT_PRIMARY = "600 16px ui-sans-serif, system-ui, sans-serif";
export const TITLE_FONT_SECONDARY = "400 12px ui-sans-serif, system-ui, sans-serif";

/** Body / legend font (also applied via Plot's `style`). */
export const BODY_FONT = "400 12px ui-sans-serif, system-ui, sans-serif";

/** Trace + peak stroke widths — wider than on-screen for printability. */
export const TRACE_STROKE_PX = 1.75;
export const PREDICTED_Q_STROKE_PX = 1.5;

/** Light palette for export (overrides dark-theme CSS vars). Values are
 *  literal CSS colors, not var(--…) — the renderer never resolves vars. */
export const LIGHT_PALETTE = {
  background:    "#ffffff",
  text:          "#1a1a1a",
  textMuted:     "#555555",
  textDim:       "#888888",
  border:        "#cccccc",
  trace:         "#1a1a1a",
  peakAuto:      "#1f5fb0",
  peakManual:    "#a0421f",
  peakExcluded:  "rgba(31, 95, 176, 0.35)",
} as const;

/**
 * Comparison palette tuned for white background. Same hue assignments as
 * `COMPARE_PALETTE` in `lib/comparison/coloring.ts` but lower OKLCH luminance
 * so contrast holds against the export's light bg. (Per spec §Risks — the
 * dark-theme palette has L = 0.76–0.80; we drop to 0.50–0.58 here.)
 *
 * If the implementer sees a hue here that fails contrast against #fff during
 * step 3 wiring, adjust THIS constant — keep on-screen → export hue mapping
 * stable.
 *
 * **Phase-offset invariant (#252).** Every entry sits >=13deg from every
 * `PHASE_PALETTE` hue, matching the on-screen `COMPARE_PALETTE` floor — since
 * PR #208 this palette also drives heatmap fill colors, so a member trace/fill
 * must never read as perceptually identical to a phase color. Entry [8] was
 * relocated 315->285 (3deg from Fd3m 318) in PR #251 r1; entries [6] 263->249
 * (1deg from Lamellar 264) and [7] 295->279 (5deg from Ia3d 300) were re-tuned
 * in #252 to mirror the on-screen palette's hue layout. `test/presets.test.ts`
 * pins the floor — keep it green if you re-tune hues here.
 */
export const COMPARE_PALETTE_LIGHT: readonly string[] = [
  "oklch(0.55 0.16  33)", // warm coral
  "oklch(0.58 0.16  77)", // gold
  "oklch(0.55 0.14 147)", // moss
  "oklch(0.54 0.13 175)", // teal
  "oklch(0.55 0.13 200)", // cyan
  "oklch(0.54 0.13 220)", // azure
  "oklch(0.55 0.12 249)", // lavender (was 263 = 1deg from Lamellar 264; #252)
  "oklch(0.50 0.14 279)", // purple (was 295 = 5deg from Ia3d 300; #252)
  "oklch(0.50 0.14 285)", // violet-purple (was 315 = 3° from Fd3m 318; #251 r1)
  "oklch(0.50 0.14 333)", // raspberry
  "oklch(0.55 0.14   3)", // pink-red
  "oklch(0.58 0.12 105)", // chartreuse
];

/** Neutral cool-gray for orphans on light bg. */
export const ORPHAN_FALLBACK_LIGHT = "oklch(0.50 0.02 270)";
