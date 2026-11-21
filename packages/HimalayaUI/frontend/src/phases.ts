export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

export const PHASE_PALETTE: Readonly<Record<KnownPhase, string>> = {
  Pn3m:      "oklch(0.520 0.100 200)",
  Im3m:      "oklch(0.480 0.140 258)",
  Ia3d:      "oklch(0.500 0.150 300)",
  Fm3m:      "oklch(0.400 0.050 230)",
  Fd3m:      "oklch(0.420 0.120 325)",
  Hexagonal: "oklch(0.530 0.140 130)",
  Lamellar:  "oklch(0.440 0.090 100)",
  Square:    "oklch(0.550 0.020 250)",
};

export const PHASE_PALETTE_DARK: Readonly<Record<KnownPhase, string>> = {
  Pn3m:      "oklch(0.740 0.100 200)",
  Im3m:      "oklch(0.720 0.120 258)",
  Ia3d:      "oklch(0.700 0.140 300)",
  Fm3m:      "oklch(0.660 0.060 230)",
  Fd3m:      "oklch(0.760 0.120 325)",
  Hexagonal: "oklch(0.760 0.140 130)",
  Lamellar:  "oklch(0.700 0.090 100)",
  Square:    "oklch(0.840 0.020 250)",
};

const FALLBACK = "oklch(0.65 0.02 270)";

export function phaseToken(phase: KnownPhase): string {
  return `--color-phase-${phase.toLowerCase()}`;
}

export function phaseColor(phase: string): string {
  const light = (PHASE_PALETTE as Record<string, string>)[phase];
  return light ? `var(${phaseToken(phase as KnownPhase)}, ${light})` : FALLBACK;
}
