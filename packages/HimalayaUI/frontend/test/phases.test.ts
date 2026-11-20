import { describe, it, expect } from "vitest";
import { phaseColor, PHASE_PALETTE, KNOWN_PHASES } from "../src/phases";
import { contrastRatio, parseOklch } from "./helpers/contrast";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const COLOR_RE = /^(#[0-9a-f]{3,8}|oklch\(.+\)|rgb\(.+\)|var\(--.+\))$/i;
const PLATE = "oklch(0.992 0.004 90)";
const MIN_DISTANCE = 0.05;
const MIN_HUE_GAP = 30;
const CHROMATIC = 0.06;

const css = readFileSync(resolve(__dirname, "../src/styles.css"), "utf8");
const design = readFileSync(resolve(__dirname, "../../../../DESIGN.md"), "utf8");

function cssToken(name: string): string {
  const m = new RegExp(`--color-${name}:\\s*(oklch\\([^)]*\\))`).exec(css);
  if (!m) throw new Error(`--color-${name} not found in styles.css`);
  return m[1]!;
}

function oklabDistance(x: string, y: string): number {
  const p = parseOklch(x);
  const q = parseOklch(y);
  const ab = (c: { C: number; h: number }) => [c.C * Math.cos((c.h * Math.PI) / 180), c.C * Math.sin((c.h * Math.PI) / 180)];
  const [pa, pb] = ab(p);
  const [qa, qb] = ab(q);
  return Math.hypot(p.L - q.L, pa! - qa!, pb! - qb!);
}

function hueGap(x: string, y: string): number {
  const d = Math.abs(parseOklch(x).h - parseOklch(y).h) % 360;
  return Math.min(d, 360 - d);
}

const STATUS = ["accent", "success", "warning", "error"].map((n) => [n, cssToken(n)] as const);
const INK = ["ink", "ink-soft", "ink-faint"].map((n) => [n, cssToken(n)] as const);

describe("phases", () => {
  it("returns a distinct color for each known phase", () => {
    const colors = new Set(KNOWN_PHASES.map((p) => phaseColor(p)));
    expect(colors.size).toBe(KNOWN_PHASES.length);
    for (const c of colors) expect(c).toMatch(COLOR_RE);
  });

  it("returns a fallback color for unknown phases", () => {
    expect(phaseColor("Unknown")).toMatch(COLOR_RE);
  });

  it("matches the phase-* values in DESIGN.md", () => {
    for (const phase of KNOWN_PHASES) {
      const m = new RegExp(`^\\s*phase-${phase.toLowerCase()}:\\s*"(oklch\\([^)]*\\))"`, "m").exec(design);
      expect(m, `phase-${phase.toLowerCase()} in DESIGN.md`).not.toBeNull();
      expect(PHASE_PALETTE[phase]).toBe(m![1]);
    }
  });

  it("every phase color passes WCAG AA (>= 4.5:1) on --plate", () => {
    for (const phase of KNOWN_PHASES) {
      const ratio = contrastRatio(PHASE_PALETTE[phase], PLATE);
      expect(ratio, `${phase} (${PHASE_PALETTE[phase]}) on --plate`).toBeGreaterThanOrEqual(4.5);
    }
  });
});

describe("phases — separation from chrome and status colours", () => {
  it("no phase colour sits within 0.05 OKLab of the accent, a status colour or an ink tone", () => {
    for (const phase of KNOWN_PHASES) {
      for (const [name, value] of [...STATUS, ...INK]) {
        expect(oklabDistance(PHASE_PALETTE[phase], value), `${phase} vs ${name}`).toBeGreaterThanOrEqual(MIN_DISTANCE);
      }
    }
  });

  it("every chromatic phase hue is at least 30° from the accent and every status hue", () => {
    for (const phase of KNOWN_PHASES) {
      if (parseOklch(PHASE_PALETTE[phase]).C < CHROMATIC) continue;
      for (const [name, value] of STATUS) {
        expect(hueGap(PHASE_PALETTE[phase], value), `${phase} vs ${name}`).toBeGreaterThanOrEqual(MIN_HUE_GAP);
      }
    }
  });

  it("every pair of phase colours is at least 0.05 OKLab apart", () => {
    for (let i = 0; i < KNOWN_PHASES.length; i++) {
      for (let j = i + 1; j < KNOWN_PHASES.length; j++) {
        const a = KNOWN_PHASES[i]!;
        const b = KNOWN_PHASES[j]!;
        expect(oklabDistance(PHASE_PALETTE[a], PHASE_PALETTE[b]), `${a} vs ${b}`).toBeGreaterThanOrEqual(MIN_DISTANCE);
      }
    }
  });
});
