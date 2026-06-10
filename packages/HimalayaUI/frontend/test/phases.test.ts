import { describe, it, expect } from "vitest";
import { phaseColor, PHASE_PALETTE, KNOWN_PHASES } from "../src/phases";
// Shared OKLCH→linear-sRGB→WCAG colour math (extracted from this file so the
// token-contrast guard reuses the exact same coefficients/formula).
import { contrastRatio } from "./helpers/contrast";

// Accepts CSS hex (#aabbcc), oklch(...), or any other valid CSS color string —
// we just want to confirm phaseColor returns a non-empty string.
const COLOR_RE = /^(#[0-9a-f]{3,8}|oklch\(.+\)|rgb\(.+\)|var\(--.+\))$/i;

// --plate from DESIGN.md (the paper plate the phase marks sit on).
const PLATE = "oklch(0.992 0.004 90)";

// DESIGN.md `colors:` block — the authoritative paper-tuned phase palette.
const DESIGN_PHASE_VALUES: Record<string, string> = {
  Pn3m: "oklch(0.570 0.150 58)",
  Im3m: "oklch(0.520 0.120 162)",
  Ia3d: "oklch(0.520 0.130 300)",
  Fm3m: "oklch(0.550 0.140 18)",
  Fd3m: "oklch(0.520 0.120 318)",
  Hexagonal: "oklch(0.540 0.130 350)",
  Lamellar: "oklch(0.505 0.150 264)",
  Square: "oklch(0.550 0.130 132)",
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
