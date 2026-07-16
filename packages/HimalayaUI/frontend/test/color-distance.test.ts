/**
 * Unit tests for the shared angular-hue-distance helper extracted from
 * coloring.test.ts (issue #252). Load-bearing for the phase-offset >=13deg
 * invariant asserted in both coloring.test.ts (on-screen COMPARE_PALETTE) and
 * presets.test.ts (export COMPARE_PALETTE_LIGHT).
 */
import { describe, it, expect } from "vitest";
import { angularHueDistance } from "../src/lib/color-distance";

describe("angularHueDistance", () => {
  it("is zero for identical hues", () => {
    expect(angularHueDistance(120, 120)).toBe(0);
  });

  it("is symmetric", () => {
    expect(angularHueDistance(10, 200)).toBe(angularHueDistance(200, 10));
  });

  it("takes the shortest arc across the 0/360 wraparound", () => {
    // 350 -> 10 is 20deg the short way, not 340deg the long way.
    expect(angularHueDistance(350, 10)).toBe(20);
    expect(angularHueDistance(10, 350)).toBe(20);
  });

  it("never exceeds 180deg", () => {
    expect(angularHueDistance(0, 270)).toBe(90);
    expect(angularHueDistance(0, 180)).toBe(180);
  });

  it("matches the known phase-offset distances (regression anchors)", () => {
    // From issue #252: light palette violations vs phase hues.
    expect(angularHueDistance(263, 264)).toBe(1);   // [6] vs Lamellar
    expect(angularHueDistance(295, 300)).toBe(5);   // [7] vs Ia3d
    // ...and the chosen replacements clear the floor (exact on-screen mirror).
    expect(angularHueDistance(249, 264)).toBe(15);  // [6]->249 vs Lamellar
    expect(angularHueDistance(279, 300)).toBe(21);  // [7]->279 vs Ia3d
  });
});
