// yBands.ts — pure y-band layout math for the on-screen multi-trace plot
// (MultiTracePlot, MemberMetaGutter). Lives in lib/ so the layout math stays
// presentation-free. (The export figure, cleanFigureSvg, computes its own band
// geometry; the former Observable Plot export mark that also consumed this was
// retired.)
//
// Originally lived in MultiTracePlot.tsx; extracted as part of issue #93's
// follow-up to close the layering smell flagged in PR #99 review.

/**
 * Compute the per-member y-band envelopes (top, bottom in pixels) given
 * the member band ratios and the panel pixel height. Pure function; the
 * order of `bandRatios` matches the caller's render order (display_order).
 *
 * Degenerate inputs:
 *   - `bandRatios` empty → `[]`.
 *   - `total ≤ 0` (every ratio zero or negative) → equal-width slicing so
 *     bands are well-defined rather than NaN.
 */
export function computeYBands(
  bandRatios: number[],
  panelHeight: number,
): Array<[number, number]> {
  if (bandRatios.length === 0) return [];
  const total = bandRatios.reduce((s, r) => s + Math.max(0, r), 0);
  if (total <= 0) {
    const each = panelHeight / bandRatios.length;
    return bandRatios.map((_, i) => [i * each, (i + 1) * each]);
  }
  const out: Array<[number, number]> = [];
  let cumulative = 0;
  for (const r of bandRatios) {
    const top    = (cumulative / total) * panelHeight;
    const next   = cumulative + Math.max(0, r);
    const bottom = (next / total) * panelHeight;
    out.push([top, bottom]);
    cumulative = next;
  }
  return out;
}

/**
 * Plate-mirroring band layout for the export waterfall (BU-EXPORTDIVERGE).
 * Mirrors WaterfallChart's stacking model, then normalizes into the fixed
 * export panel:
 *
 *   - display-order row 0 sits at the BOTTOM (the plate paints bottom-up;
 *     `computeYBands` above paints top-down — using it published a vertically
 *     flipped figure);
 *   - the inter-row STEP scales with `offsetScale` (the builder's offset
 *     slider, same 0.1 floor as the plate) while each row's band height stays
 *     constant, so separation/overlap ratios match the on-screen figure;
 *   - per-member durable `y_offset` nudges are honored, like the plate;
 *   - the resulting stack is rescaled to exactly fill `panelHeight`, because
 *     the export canvas is fixed while the on-screen stack grows with the
 *     offset. Relative geometry — the WYSIWYG axis — is what's preserved.
 *
 * Returned bands are indexed by member position (same convention as
 * `computeYBands`), each `[top, bottom]` in panel pixels.
 */
export function computeWaterfallExportBands(
  rows: ReadonlyArray<{ band_height: number; y_offset: number }>,
  panelHeight: number,
  offsetScale: number,
): Array<[number, number]> {
  const n = rows.length;
  if (n === 0) return [];
  const scale = Math.max(0.1, offsetScale); // plate's collapse floor
  const totalWeight = rows.reduce((s, r) => s + Math.max(0, r.band_height), 0) || n;
  // y_offset is a durable nudge in PLATE pixels (WaterfallChart's 420px
  // internal stack). Band heights are weight-normalized in each space, but a
  // raw-pixel nudge would carry ~1.5x more relative weight in the smaller
  // export panel - convert it to panel space so the nudge:band ratio matches
  // the plate.
  const PLATE_TOTAL_H = 420;
  const nudgeScale = panelHeight / PLATE_TOTAL_H;
  const bands: Array<[number, number]> = new Array<[number, number]>(n);
  let cumulative = 0;
  let stackHeight = 0;
  // Reversed walk = top-down pixel placement of a bottom-up display order,
  // exactly as WaterfallChart places its rows.
  for (let i = n - 1; i >= 0; i--) {
    const r = rows[i]!;
    const h = (Math.max(0, r.band_height) / totalWeight) * panelHeight;
    const top = cumulative + r.y_offset * nudgeScale;
    cumulative += h * scale;
    stackHeight = Math.max(stackHeight, top + h);
    bands[i] = [top, top + h];
  }
  const k = panelHeight / Math.max(stackHeight, 1e-9);
  return bands.map(([top, bottom]) => [top * k, bottom * k]);
}
