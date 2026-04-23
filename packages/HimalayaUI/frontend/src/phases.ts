export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

const PALETTE: Record<KnownPhase, string> = {
  Pn3m:      "#d97706",
  Im3m:      "#22c55e",
  Ia3d:      "#3b82f6",
  Fm3m:      "#a855f7",
  Fd3m:      "#ec4899",
  Hexagonal: "#eab308",
  Lamellar:  "#14b8a6",
  Square:    "#f97316",
};

const FALLBACK = "#9a9894";

export function phaseColor(phase: string): string {
  return (PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
