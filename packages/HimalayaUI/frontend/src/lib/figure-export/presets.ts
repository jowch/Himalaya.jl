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
export const PEAK_TICK_STROKE_PX = 1.5;
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
 * **Phase-offset parity (PR #251 r1).** Entry [8] was relocated from hue 315
 * to 285 to keep parity with `COMPARE_PALETTE`'s same-index fix (3° collision
 * with Fd3m 318). Several other entries in this light variant — 175, 263,
 * 295, 3 — sit closer to phase hues than the on-screen ≥13° floor; tracked
 * but out of scope for #208 (issue covers the on-screen render-core).
 */
export const COMPARE_PALETTE_LIGHT: readonly string[] = [
  "oklch(0.55 0.16  33)", // warm coral
  "oklch(0.58 0.16  77)", // gold
  "oklch(0.55 0.14 147)", // moss
  "oklch(0.54 0.13 175)", // teal
  "oklch(0.55 0.13 200)", // cyan
  "oklch(0.54 0.13 220)", // azure
  "oklch(0.55 0.12 263)", // lavender
  "oklch(0.50 0.14 295)", // purple
  "oklch(0.50 0.14 285)", // violet-purple (was 315 = 3° from Fd3m 318; #251 r1)
  "oklch(0.50 0.14 333)", // raspberry
  "oklch(0.55 0.14   3)", // pink-red
  "oklch(0.58 0.12 105)", // chartreuse
];

/** Neutral cool-gray for orphans on light bg. */
export const ORPHAN_FALLBACK_LIGHT = "oklch(0.50 0.02 270)";
