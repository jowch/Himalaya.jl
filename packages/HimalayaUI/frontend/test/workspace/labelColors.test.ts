import { describe, expect, it } from "vitest";
import { LABEL_COLORS } from "../../src/api";
import { contrastRatio, oklchToLinearSrgb, parseOklch } from "../helpers/contrast";
import { themeToken } from "../helpers/cssBlocks";

const DOT_FLOOR = 0.05;
const WASH_FLOOR = 0.025;
const GAMUT_SLACK = 1e-3;

function oklab(s: string): [number, number, number] {
  const { L, C, h } = parseOklch(s);
  return [L, C * Math.cos((h * Math.PI) / 180), C * Math.sin((h * Math.PI) / 180)];
}

function oklabDistance(x: string, y: string): number {
  const [pl, pa, pb] = oklab(x);
  const [ql, qa, qb] = oklab(y);
  return Math.hypot(pl - ql, pa - qa, pb - qb);
}

function dotDistance(x: string, y: string): number {
  const [pl, pa, pb] = oklab(x);
  const [ql, qa, qb] = oklab(y);
  return Math.hypot(pl - ql, (pa - qa) / 2, (pb - qb) / 2);
}

function pairs<T>(xs: readonly T[]): [T, T][] {
  return xs.flatMap((a, i) => xs.slice(i + 1).map((b) => [a, b] as [T, T]));
}

describe.each(["light", "dark"] as const)("label colours in the %s theme", (theme) => {
  const ink = (c: string): string => themeToken(theme, `label-${c}-ink`);
  const wash = (c: string): string => themeToken(theme, `label-${c}-wash`);

  it("states an ink and a wash for each of the eight colours", () => {
    expect(LABEL_COLORS).toHaveLength(8);
    for (const c of LABEL_COLORS) expect([ink(c), wash(c)].every((v) => v.startsWith("oklch("))).toBe(true);
  });

  it("keeps every ink readable on its own wash and on the plate", () => {
    const plate = themeToken(theme, "plate");
    for (const c of LABEL_COLORS) {
      expect(contrastRatio(ink(c), wash(c)), `${c} on its wash`).toBeGreaterThanOrEqual(4.5);
      expect(contrastRatio(ink(c), plate), `${c} on plate`).toBeGreaterThanOrEqual(4.5);
    }
  });

  it("keeps every ink and wash inside sRGB, so the colour on screen is the one measured here", () => {
    for (const c of LABEL_COLORS) {
      for (const v of [ink(c), wash(c)]) {
        for (const channel of oklchToLinearSrgb(parseOklch(v))) {
          expect(channel, `${c} ${v}`).toBeGreaterThanOrEqual(-GAMUT_SLACK);
          expect(channel, `${c} ${v}`).toBeLessThanOrEqual(1 + GAMUT_SLACK);
        }
      }
    }
  });

  it("keeps every pair of inks apart as small dots, where lightness counts twice as much as hue", () => {
    for (const [a, b] of pairs(LABEL_COLORS)) expect(dotDistance(ink(a), ink(b)), `${a} vs ${b}`).toBeGreaterThanOrEqual(DOT_FLOOR);
  });

  it("keeps every pair of washes apart, so chips differ without reading their text", () => {
    for (const [a, b] of pairs(LABEL_COLORS)) expect(oklabDistance(wash(a), wash(b)), `${a} vs ${b}`).toBeGreaterThanOrEqual(WASH_FLOOR);
  });
});
