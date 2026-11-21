import { describe, expect, it } from "vitest";
import { LABEL_COLORS } from "../../src/api";
import { KNOWN_PHASES, phaseToken } from "../../src/phases";
import { contrastRatio, parseOklch } from "../helpers/contrast";
import { themeToken } from "../helpers/cssBlocks";

const FROM_OTHERS = 0.06;
const BETWEEN_LABELS = 0.047;

function oklabDistance(x: string, y: string): number {
  const ab = (s: string): [number, number, number] => {
    const { L, C, h } = parseOklch(s);
    return [L, C * Math.cos((h * Math.PI) / 180), C * Math.sin((h * Math.PI) / 180)];
  };
  const [pl, pa, pb] = ab(x);
  const [ql, qa, qb] = ab(y);
  return Math.hypot(pl - ql, pa - qa, pb - qb);
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

  it("keeps every ink at least 0.06 OKLab from each phase colour and status colour", () => {
    const others = [
      ...KNOWN_PHASES.map((p) => [p, themeToken(theme, phaseToken(p).replace("--color-", ""))] as const),
      ...["accent", "success", "warning", "error"].map((n) => [n, themeToken(theme, n)] as const),
    ];
    for (const c of LABEL_COLORS) {
      for (const [name, value] of others) expect(oklabDistance(ink(c), value), `${c} vs ${name}`).toBeGreaterThanOrEqual(FROM_OTHERS);
    }
  });

  it("keeps every pair of inks at least 0.047 OKLab apart", () => {
    for (let i = 0; i < LABEL_COLORS.length; i++) {
      for (let j = i + 1; j < LABEL_COLORS.length; j++) {
        const a = LABEL_COLORS[i]!;
        const b = LABEL_COLORS[j]!;
        expect(oklabDistance(ink(a), ink(b)), `${a} vs ${b}`).toBeGreaterThanOrEqual(BETWEEN_LABELS);
      }
    }
  });
});
