/**
 * Shortest angular distance between two hues, in degrees (mod 360, range
 * [0, 180]). Extracted from `test/coloring.test.ts` (issue #252) so the
 * phase-offset >=13deg invariant can run against both the on-screen
 * `COMPARE_PALETTE` (`lib/comparison/coloring.ts`) and the export
 * `COMPARE_PALETTE_LIGHT` (`lib/figure-export/presets.ts`) from one proven
 * implementation. Phase/member colors visually conflate when their hues sit
 * too close; this is the load-bearing perceptual check behind that floor.
 */
export function angularHueDistance(a: number, b: number): number {
  const d = (Math.abs(((a - b) % 360) + 540) % 360) - 180;
  return Math.abs(d);
}
