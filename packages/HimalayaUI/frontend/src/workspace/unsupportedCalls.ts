import type { IndexEntry, SampleCallMember } from "../api";

export function unsupportedCalls(
  members: readonly SampleCallMember[],
  offered: readonly IndexEntry[],
  lattices: readonly IndexEntry[],
): IndexEntry[] {
  const phases = new Set(offered.map((c) => c.phase));
  const byId = new Map(lattices.map((ix) => [ix.id, ix]));
  const out: IndexEntry[] = [];
  for (const m of members) {
    if (phases.has(m.phase)) continue;
    phases.add(m.phase);
    out.push(byId.get(m.index_id) ?? placeholder(m));
  }
  return out;
}

function placeholder(m: SampleCallMember): IndexEntry {
  return {
    id: m.index_id, exposure_id: m.exposure_id, phase: m.phase, basis: 0, score: null, r_squared: null,
    lattice_d: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}
