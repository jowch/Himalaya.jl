// curvature.ts — phase-averaged Gaussian curvature κ for the three bicontinuous
// cubic phases, computed client-side from a fitted lattice parameter. This is a
// faithful port of the backend `_ngc_for_phase` (packages/HimalayaUI/src/
// routes_analysis.jl): the snapshot already carries κ on the *confirmed_index*
// (`ngc`), but a SECOND coexisting phase (MemberSnapshotPhase) carries only its
// lattice_d, so the figure recomputes κ from (phase, lattice_d) for any phase —
// same formula, same constants, so the two never disagree.
//
// Formula: κ = -2π·(χ/A₀)/a². Units are 1/length² in the lattice's own unit
// (Å⁻² when lattice_d is in Å, nm⁻² when in nm) — no hard-coded unit conversion,
// matching the backend's unit-naked form. κ is defined ONLY for the three
// bicontinuous cubics (Pn3m / Im3m / Ia3d); every other phase returns null.

/** Euler characteristic χ and dimensionless surface area A₀ per unit cell, for
 *  the three bicontinuous cubic minimal surfaces (index.jl `ngc` constants). */
const CUBIC_CURVATURE: Record<string, { chi: number; A0: number }> = {
  Pn3m: { chi: -2, A0: 1.919 },
  Im3m: { chi: -4, A0: 2.345 },
  Ia3d: { chi: -8, A0: 3.091 },
};

/** Phase-averaged Gaussian curvature κ, or null when the phase is not a
 *  bicontinuous cubic or the lattice is missing / non-positive. */
export function kappaForPhase(phase: string, latticeD: number | null): number | null {
  if (latticeD == null || !(latticeD > 0)) return null;
  const c = CUBIC_CURVATURE[phase];
  if (!c) return null;
  return (-2 * Math.PI * (c.chi / c.A0)) / (latticeD * latticeD);
}
