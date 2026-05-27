import { describe, it, expect } from "vitest";
import { phaseColor, PHASE_PALETTE, KNOWN_PHASES } from "../src/phases";

// Accepts CSS hex (#aabbcc), oklch(...), or any other valid CSS color string —
// we just want to confirm phaseColor returns a non-empty string.
const COLOR_RE = /^(#[0-9a-f]{3,8}|oklch\(.+\)|rgb\(.+\)|var\(--.+\))$/i;

// --plate from DESIGN.md (the paper plate the phase marks sit on).
const PLATE = "oklch(0.992 0.004 90)";

// --- Minimal OKLCH → linear sRGB → WCAG relative luminance / contrast ---------
// Self-contained so the test pulls in no color library. Implements the OKLab
// inverse transform (Björn Ottosson) + the standard sRGB OETF + WCAG 2.x
// relative-luminance + contrast-ratio formulas.

function parseOklch(s: string): { L: number; C: number; h: number } {
  const m = /oklch\(\s*([\d.]+)\s+([\d.]+)\s+([\d.-]+)\s*\)/i.exec(s);
  if (!m) throw new Error(`not an oklch() string: ${s}`);
  return { L: parseFloat(m[1]!), C: parseFloat(m[2]!), h: parseFloat(m[3]!) };
}

function oklchToLinearSrgb({ L, C, h }: { L: number; C: number; h: number }): [number, number, number] {
  const hr = (h * Math.PI) / 180;
  const a = C * Math.cos(hr);
  const b = C * Math.sin(hr);

  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.291485548 * b;

  const l = l_ * l_ * l_;
  const m = m_ * m_ * m_;
  const s = s_ * s_ * s_;

  return [
    +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s,
  ];
}

// WCAG relative luminance expects *linear* channel values (which is exactly
// what the OKLab inverse produces before applying the sRGB OETF), clamped to
// [0,1] for out-of-gamut colors.
function relativeLuminance(s: string): number {
  const [r, g, b] = oklchToLinearSrgb(parseOklch(s)).map((c) => Math.min(1, Math.max(0, c)));
  return 0.2126 * r! + 0.7152 * g! + 0.0722 * b!;
}

function contrastRatio(a: string, b: string): number {
  const la = relativeLuminance(a);
  const lb = relativeLuminance(b);
  const lighter = Math.max(la, lb);
  const darker = Math.min(la, lb);
  return (lighter + 0.05) / (darker + 0.05);
}

// DESIGN.md `colors:` block — the authoritative paper-tuned phase palette.
const DESIGN_PHASE_VALUES: Record<string, string> = {
  Pn3m: "oklch(0.585 0.150 58)",
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
