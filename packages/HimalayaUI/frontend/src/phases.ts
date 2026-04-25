export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

// Earthy / muted OKLCH palette. Indices are *context*, not the active edit —
// they should sit quietly behind the trace and the peak markers. Chroma is
// kept low (0.07–0.10) and luminance roughly equal across phases so no single
// hue dominates by brightness. Hues are spaced to stay distinguishable while
// avoiding ±25° of either peak color (--color-accent at 220°,
// --color-peak-manual at 340°) and the warning yellow (~75° at high chroma).
const PALETTE: Record<KnownPhase, string> = {
  Pn3m:      "oklch(0.74 0.10  35)", // terracotta
  Im3m:      "oklch(0.74 0.08 145)", // sage
  Ia3d:      "oklch(0.70 0.09 280)", // periwinkle (clear of 220° accent)
  Fm3m:      "oklch(0.72 0.09 320)", // mauve
  Fd3m:      "oklch(0.72 0.10   5)", // dusty rose
  Hexagonal: "oklch(0.78 0.09  75)", // ochre (low-chroma, not warning yellow)
  Lamellar:  "oklch(0.74 0.07 180)", // dusty teal
  Square:    "oklch(0.74 0.10  55)", // amber-tan
};

const FALLBACK = "oklch(0.65 0.02 270)"; // neutral cool-gray

export function phaseColor(phase: string): string {
  return (PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
