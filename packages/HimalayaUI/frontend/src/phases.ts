export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

export const CUBIC_PHASES: ReadonlySet<string> = new Set([
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
]);

export const PHASE_PALETTE: Readonly<Record<KnownPhase, string>> = {
  Pn3m:      "oklch(0.520 0.100 200)",
  Im3m:      "oklch(0.480 0.140 258)",
  Ia3d:      "oklch(0.500 0.150 300)",
  Fm3m:      "oklch(0.400 0.050 230)",
  Fd3m:      "oklch(0.420 0.120 325)",
  Hexagonal: "oklch(0.500 0.130 138)",
  Lamellar:  "oklch(0.500 0.100 105)",
  Square:    "oklch(0.480 0.030 80)",
};

const FALLBACK = "oklch(0.65 0.02 270)"; // neutral cool-gray

export function phaseColor(phase: string): string {
  return (PHASE_PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
