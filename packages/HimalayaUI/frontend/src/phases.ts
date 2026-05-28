export const KNOWN_PHASES = [
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
  "Hexagonal", "Lamellar", "Square",
] as const;

export type KnownPhase = typeof KNOWN_PHASES[number];

export const CUBIC_PHASES: ReadonlySet<string> = new Set([
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
]);

// Paper-tuned OKLCH phase palette — "The Print" (R0b, finding A-2).
// Authoritative source: DESIGN.md `colors:` block (phase-* values). These are
// the exact values used here; the DESIGN.md-match test in test/phases.test.ts
// pins them, so update both together.
//
// Tuning: luminance L≈0.50–0.58, chroma 0.12–0.15. Darker + higher-chroma than
// the retired dark-field palette (L 0.76–0.80) so each hue clears WCAG AA
// (≥4.5:1) on the light --plate paper (oklch(0.992 0.004 90)); test/phases.test.ts
// asserts the AA floor per phase.
//
// Colour-blind second channel: every phase hue is *always* paired with the
// phase name/label (and shape, e.g. the Miller-plot symbols) in the consuming
// components (PhasePanel chip, MentionChip, MemberTraceLayer) — colour is never
// the sole signal. This palette only supplies the hue; do not make any consumer
// colour-only.
//
// Exported so the Compare page's categorical trace palette
// (`COMPARE_PALETTE` in `lib/comparison/coloring.ts`) can assert it does
// not collide — `byPhase` and `bySample` would visually conflate if they
// shared colors.
export const PHASE_PALETTE: Readonly<Record<KnownPhase, string>> = {
  Pn3m:      "oklch(0.570 0.150 58)",  // amber — L nudged 0.585→0.570 to clear AA on --plate (see note)
  Im3m:      "oklch(0.520 0.120 162)", // sage
  Ia3d:      "oklch(0.520 0.130 300)", // violet
  Fm3m:      "oklch(0.550 0.140 18)",  // coral
  Fd3m:      "oklch(0.520 0.120 318)", // rose-purple
  Hexagonal: "oklch(0.540 0.130 350)", // magenta-rose
  Lamellar:  "oklch(0.505 0.150 264)", // indigo
  Square:    "oklch(0.550 0.130 132)", // green
};

const FALLBACK = "oklch(0.65 0.02 270)"; // neutral cool-gray

export function phaseColor(phase: string): string {
  return (PHASE_PALETTE as Record<string, string>)[phase] ?? FALLBACK;
}
