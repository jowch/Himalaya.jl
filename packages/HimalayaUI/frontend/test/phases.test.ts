import { describe, it, expect } from "vitest";
import { phaseColor, PHASE_PALETTE, KNOWN_PHASES } from "../src/phases";
// Shared OKLCH→linear-sRGB→WCAG colour math (extracted from this file so the
// token-contrast guard reuses the exact same coefficients/formula).
import { contrastRatio, parseOklch } from "./helpers/contrast";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

// Accepts CSS hex (#aabbcc), oklch(...), or any other valid CSS color string —
// we just want to confirm phaseColor returns a non-empty string.
const COLOR_RE = /^(#[0-9a-f]{3,8}|oklch\(.+\)|rgb\(.+\)|var\(--.+\))$/i;

// --plate from DESIGN.md (the paper plate the phase marks sit on).
const PLATE = "oklch(0.992 0.004 90)";

// DESIGN.md `colors:` block — the authoritative paper-tuned phase palette.
const DESIGN_PHASE_VALUES: Record<string, string> = {
  Pn3m: "oklch(0.520 0.100 200)",
  Im3m: "oklch(0.480 0.140 258)",
  Ia3d: "oklch(0.500 0.150 300)",
  Fm3m: "oklch(0.400 0.050 230)",
  Fd3m: "oklch(0.420 0.120 325)",
  Hexagonal: "oklch(0.500 0.130 138)",
  Lamellar: "oklch(0.500 0.100 105)",
  Square: "oklch(0.480 0.030 80)",
};

describe("phases", () => {
  it("returns a distinct color for each known phase", () => {
    const colors = new Set(KNOWN_PHASES.map((p) => phaseColor(p)));
    expect(colors.size).toBe(KNOWN_PHASES.length);
    for (const c of colors) expect(c).toMatch(COLOR_RE);
  });

  it("returns a fallback color for unknown phases", () => {
    expect(phaseColor("Unknown")).toMatch(COLOR_RE);
  });
});

describe("phases — paper-tuned palette (R0b / finding A-2)", () => {
  it("matches the DESIGN.md authoritative phase values exactly", () => {
    for (const phase of KNOWN_PHASES) {
      expect(PHASE_PALETTE[phase]).toBe(DESIGN_PHASE_VALUES[phase]);
    }
  });

  it("every phase color passes WCAG AA (>= 4.5:1) on --plate", () => {
    for (const phase of KNOWN_PHASES) {
      const ratio = contrastRatio(PHASE_PALETTE[phase], PLATE);
      expect(ratio, `${phase} (${PHASE_PALETTE[phase]}) on --plate`).toBeGreaterThanOrEqual(4.5);
    }
  });
});

function oklabDistance(x: string, y: string): number {
  const p = parseOklch(x);
  const q = parseOklch(y);
  const ab = (c: { C: number; h: number }) => [c.C * Math.cos((c.h * Math.PI) / 180), c.C * Math.sin((c.h * Math.PI) / 180)];
  const [pa, pb] = ab(p);
  const [qa, qb] = ab(q);
  return Math.hypot(p.L - q.L, pa! - qa!, pb! - qb!);
}

describe("phases — separation from chrome and status colours", () => {
  const css = readFileSync(resolve(__dirname, "../src/styles.css"), "utf8");
  const tokens = ["accent", "success", "warning", "error"].map((name) => {
    const m = new RegExp(`--color-${name}:\\s*(oklch\\([^)]*\\))`).exec(css);
    if (!m) throw new Error(`--color-${name} not found in styles.css`);
    return [name, m[1]!] as const;
  });

  it("reads every status token from styles.css", () => {
    expect(tokens).toHaveLength(4);
  });

  it("no phase colour sits within 0.05 OKLab of the accent or a status colour", () => {
    for (const phase of KNOWN_PHASES) {
      for (const [name, value] of tokens) {
        const d = oklabDistance(PHASE_PALETTE[phase], value);
        expect(d, `${phase} vs ${name}`).toBeGreaterThanOrEqual(0.05);
      }
    }
  });

  it("every pair of phase colours is at least 0.05 OKLab apart", () => {
    for (let i = 0; i < KNOWN_PHASES.length; i++) {
      for (let j = i + 1; j < KNOWN_PHASES.length; j++) {
        const a = KNOWN_PHASES[i]!;
        const b = KNOWN_PHASES[j]!;
        expect(oklabDistance(PHASE_PALETTE[a], PHASE_PALETTE[b]), `${a} vs ${b}`).toBeGreaterThanOrEqual(0.05);
      }
    }
  });
});
