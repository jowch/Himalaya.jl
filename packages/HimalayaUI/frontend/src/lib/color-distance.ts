/**
 * Shortest angular distance between two hues, in degrees (mod 360, range
 * [0, 180]). Extracted from `test/coloring.test.ts` (issue #252) so the
 * phase-offset >=13deg invariant can run against the on-screen
 * `COMPARE_PALETTE` (`lib/comparison/coloring.ts`) from one proven
 * implementation. Phase/member colors visually conflate when their hues sit
 * too close; this is the load-bearing perceptual check behind that floor.
 *
 * Body is the verbatim test-local helper from `coloring.test.ts`; the added
 * parens around `Math.abs(...) % 360` are a readability clarification only —
 * `%` already binds tighter than `-`, so the grouping was implicit and the
 * output is identical (pinned by `test/color-distance.test.ts` anchors).
 */
export function angularHueDistance(a: number, b: number): number {
  const d = (Math.abs(((a - b) % 360) + 540) % 360) - 180;
  return Math.abs(d);
}
