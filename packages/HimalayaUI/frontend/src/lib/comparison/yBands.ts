// yBands.ts — pure y-band layout math shared between the on-screen plot
// (MultiTracePlot, MemberMetaGutter) and the export mark builder
// (lib/figure-export/marks/multiTraceExportMarks). Lives in lib/ so the
// figure-export module doesn't have to reach into components/.
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
