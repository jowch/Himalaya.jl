/**
 * Figure-export palette invariant (issue #252).
 *
 * `COMPARE_PALETTE_LIGHT` is the light-background sibling of the on-screen
 * `COMPARE_PALETTE`. Since PR #208 it also drives heatmap fill colors in
 * figure export, so a member trace/fill must never read as perceptually
 * identical to a phase color — the same >=13deg hue-offset floor the on-screen
 * palette holds (see test/coloring.test.ts) applies here.
 *
 * Floor rationale: at L~0.55 C~0.13 a 13deg hue shift is dE2000 ~ 2.5 — the
 * smallest gap that keeps byPhase vs bySample/distinct from conflating. The
 * eight phase hues pack the warm + purple sectors tightly enough that >=15deg
 * everywhere AND twelve mutually-distinct entries is infeasible, so 13deg (not
 * 15deg) is the floor. Keep this floor numeric and in lockstep with both
 * palettes' docstrings if hues are re-laid.
 */
import { describe, it, expect } from "vitest";
import { COMPARE_PALETTE_LIGHT } from "../src/lib/figure-export/presets";
import { PHASE_PALETTE } from "../src/phases";
import { angularHueDistance } from "../src/lib/color-distance";

/** Parse the hue channel out of an `oklch(L C h)` string. */
function oklchHue(s: string): number {
  const m = /oklch\(\s*[\d.]+\s+[\d.]+\s+([\d.-]+)\s*\)/i.exec(s);
  if (!m) throw new Error(`not an oklch() string: ${s}`);
  return parseFloat(m[1]!);
}

const PHASE_OFFSET_FLOOR_DEG = 13;

describe("COMPARE_PALETTE_LIGHT — phase-offset invariant (#252)", () => {
  it("exports a palette of 10-12 colors", () => {
    expect(COMPARE_PALETTE_LIGHT.length).toBeGreaterThanOrEqual(10);
    expect(COMPARE_PALETTE_LIGHT.length).toBeLessThanOrEqual(12);
  });

  it("no entry exactly equals a phase color (string-equality first-line check)", () => {
    const phaseColors = new Set(Object.values(PHASE_PALETTE));
    for (const c of COMPARE_PALETTE_LIGHT) {
      expect(phaseColors.has(c), `${c} collides with a PHASE_PALETTE entry`).toBe(false);
    }
  });

  it("every entry sits >=13deg from every PHASE_PALETTE hue", () => {
    const phaseHues = Object.values(PHASE_PALETTE).map(oklchHue);
    for (const c of COMPARE_PALETTE_LIGHT) {
      const ch = oklchHue(c);
      for (const ph of phaseHues) {
        const dist = angularHueDistance(ch, ph);
        expect(
          dist,
          `LIGHT ${c} (hue ${ch}) vs PHASE hue ${ph}`,
        ).toBeGreaterThanOrEqual(PHASE_OFFSET_FLOOR_DEG);
      }
    }
  });
});
