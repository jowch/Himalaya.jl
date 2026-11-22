// physics: Gauss–Bonnet over one unit cell of the underlying minimal surface,
// physics: ∫K dA = 2πχ, with cell area A₀a², gives ⟨K⟩ = 2πχ / (A₀a²) (Hyde et al., The
// physics: Language of Shape, 1997). χ and A₀ are per cubic cell of the indexed lattice:
// physics: P (Im3m) −4, 2.345 and G (Ia3d) −8, 3.091 on the conventional cell; D (Pn3m)
// physics: −2, 1.919 on the primitive Pn3m cell, which is the conventional diamond cell's
// physics: −16, 3.838 at half its edge. Their √(χ/A₀) ratios reproduce the
// physics: 1 : 1.279 : 1.576 Bonnet lattice ratios in src/bonnet.jl.
const SURFACES: Readonly<Record<string, { chi: number; area: number }>> = {
  Pn3m: { chi: -2, area: 1.919 },
  Im3m: { chi: -4, area: 2.345 },
  Ia3d: { chi: -8, area: 3.091 },
};

export function gaussianCurvature(phase: string, latticeAngstrom: number | null): number | null {
  const s = SURFACES[phase];
  if (s === undefined || latticeAngstrom === null || !(latticeAngstrom > 0)) return null;
  const a = latticeAngstrom / 10;
  return (2 * Math.PI * s.chi) / (s.area * a * a);
}

export function formatCurvature(k: number): string {
  return `κ ${k.toFixed(4).replace("-", "−")} nm⁻²`;
}
