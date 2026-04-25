export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

// Vivid OKLCH palette inspired by the v3 design reference (phaseA sage/phaseB amber/phaseC violet).
// Chroma ~0.12–0.13, luminance 0.76–0.80 for clear colour on dark backgrounds.
// Hues avoid ±20° of accent (220°) and peak-manual (340°) and the high-chroma warning zone (~75°).
// Exception: Pn3m at 62° is 13° from the warning hue, but is visually distinct because its
// luminance (0.80) and chroma (0.13) differ enough from --color-warning (0.76 / 0.11 / 75°).
const PALETTE: Record<KnownPhase, string> = {
  Pn3m:      "oklch(0.80 0.13  62)", // amber — 13° from warning zone; distinct by luminance + chroma
  Im3m:      "oklch(0.78 0.12 160)", // sage   (mock-up phaseA)
  Ia3d:      "oklch(0.76 0.12 300)", // violet (mock-up phaseC)
  Fm3m:      "oklch(0.78 0.13  18)", // coral
  Fd3m:      "oklch(0.76 0.12 318)", // rose-purple (clear of 340° manual-peak)
  Hexagonal: "oklch(0.79 0.12 185)", // seafoam teal
  Lamellar:  "oklch(0.80 0.10 248)", // periwinkle (28° from 220° accent)
  Square:    "oklch(0.79 0.12 132)", // chartreuse
};

const FALLBACK = "oklch(0.65 0.02 270)"; // neutral cool-gray

export function phaseColor(phase: string): string {
  return (PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
