/**
 * series/seriesReading — derive the "phases present" reading (Plan E, Task E-4).
 *
 * The Series rail's job is to answer "what changes across the variable?" in
 * words, derived purely from per-member snapshots — NO hand-written "X → Y"
 * narration (which wouldn't generalize). For each phase that appears in any
 * member's `confirmed_phases`, we report the variable span it covers and its
 * lattice trend (first → last fitted lattice, collapsed to a single value when
 * constant). We also surface the variable values where >1 phase coexists, and
 * those that are form-factor-only.
 *
 * Cubic phases report `a` (lattice parameter); lamellar/hexagonal report `d`.
 * A `null`-state member (featureless, deliberately un-assigned) contributes
 * nothing — it is neither an indexed phase nor a form-factor reading.
 */
import type { SeriesMember, MemberSnapshotPhase } from "../../api";
import { CUBIC_PHASES } from "../../phases";

export interface SeriesReading {
  phases: Array<{ phase: string; spanLabel: string; latticeTrend: string }>;
  /** Variable values with >1 assigned phase (coexistence). */
  coexistenceAt: string[];
  /** Variable values that are form-factor-only (broad shoulder, no Bragg). */
  formFactorOnlyAt: string[];
}

/** Lattice symbol: cubics report the lattice parameter `a`; everything else
 *  (lamellar, hexagonal) reports the d-spacing `d`. */
function latticeSymbol(phase: string): "a" | "d" {
  return CUBIC_PHASES.has(phase) ? "a" : "d";
}

function fmtLattice(v: number | null): string {
  return v == null ? "—" : (Number.isInteger(v) ? String(v) : v.toFixed(1));
}

/** Per-member assignment phases, honouring the 3-state. `null` members carry no
 *  phases; form-factor members carry no phases either (but DO surface in the
 *  form-factor line). Falls back to `confirmed_index` when `confirmed_phases`
 *  is absent (older snapshots). */
function phasesOf(m: SeriesMember): MemberSnapshotPhase[] {
  const snap = m.snapshot;
  if (!snap) return [];
  if (snap.assignment_state === "null" || snap.assignment_state === "form_factor") {
    return [];
  }
  if (snap.confirmed_phases && snap.confirmed_phases.length > 0) {
    return snap.confirmed_phases;
  }
  if (snap.confirmed_index) {
    return [{ phase: snap.confirmed_index.phase, lattice_d: snap.confirmed_index.lattice_d }];
  }
  return [];
}

export function seriesReading(
  members: SeriesMember[],
  variableOf: (m: SeriesMember) => number | string,
): SeriesReading {
  // phase → ordered list of { variable, lattice } across carrying members.
  const byPhase = new Map<string, Array<{ v: string; lat: number | null }>>();
  const coexistenceAt: string[] = [];
  const formFactorOnlyAt: string[] = [];

  for (const m of members) {
    const v = String(variableOf(m));
    const phases = phasesOf(m);
    const state = m.snapshot?.assignment_state;

    if (phases.length === 0) {
      // Form-factor → surface the variable; null → contribute nothing.
      if (state === "form_factor") formFactorOnlyAt.push(v);
      continue;
    }
    if (phases.length > 1) coexistenceAt.push(v);
    for (const p of phases) {
      if (!byPhase.has(p.phase)) byPhase.set(p.phase, []);
      byPhase.get(p.phase)!.push({ v, lat: p.lattice_d });
    }
  }

  const phases = Array.from(byPhase.entries()).map(([phase, items]) => {
    const v0 = items[0]!.v;
    const v1 = items[items.length - 1]!.v;
    const l0 = items[0]!.lat;
    const l1 = items[items.length - 1]!.lat;
    const sym = latticeSymbol(phase);
    const spanLabel = v0 === v1 ? v0 : `${v0} → ${v1}`;
    const latticeTrend = l0 === l1
      ? `${sym} = ${fmtLattice(l0)} Å`
      : `${sym} ${fmtLattice(l0)} → ${fmtLattice(l1)} Å`;
    return { phase, spanLabel, latticeTrend };
  });

  return { phases, coexistenceAt, formFactorOnlyAt };
}
