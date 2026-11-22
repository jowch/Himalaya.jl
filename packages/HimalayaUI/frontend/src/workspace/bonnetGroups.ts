import type { IndexEntry } from "../api";

// physics: Bonnet-related D/P/G surfaces share Gaussian curvature κ = 2πχ/(A₀a²), so equal κ fixes a_Pn3m : a_Im3m : a_Ia3d = 1 : 1.279 : 1.576.
const BONNET_SCALE: Readonly<Record<string, number>> = { Pn3m: 1, Im3m: 1.279, Ia3d: 1.576 };
const RELTOL = 0.01;
// physics: a κ match alone admits lattice aliases, so every member must explain peaks no other member does.
const MIN_OWN = 2;

export type CallRow =
  | { kind: "one"; candidate: IndexEntry }
  | { kind: "group"; members: readonly IndexEntry[]; peaks: number };

const peakKey = (q: number) => q.toFixed(6);
const peakSet = (ix: IndexEntry) => new Set(ix.peaks.map((p) => peakKey(p.q_observed)));

function related(a: IndexEntry, b: IndexEntry): boolean {
  const sa = BONNET_SCALE[a.phase];
  const sb = BONNET_SCALE[b.phase];
  if (sa === undefined || sb === undefined || a.phase === b.phase || a.lattice_d === null || b.lattice_d === null) return false;
  const predicted = a.lattice_d * (sb / sa);
  return Math.abs(b.lattice_d - predicted) <= RELTOL * predicted;
}

function union(sets: readonly Set<string>[]): Set<string> {
  const out = new Set<string>();
  for (const s of sets) for (const k of s) out.add(k);
  return out;
}

export function bonnetGroups(candidates: readonly IndexEntry[]): IndexEntry[][] {
  const cubic = candidates
    .filter((c) => BONNET_SCALE[c.phase] !== undefined && c.lattice_d !== null && c.lattice_d > 0)
    .sort((a, b) => a.id - b.id);
  const peaks = new Map(cubic.map((c) => [c.id, peakSet(c)]));
  const of = (c: IndexEntry) => peaks.get(c.id)!;

  const groups: IndexEntry[][] = [];
  cubic.forEach((a, i) => {
    for (const b of cubic.slice(i + 1)) {
      if (!related(a, b)) continue;
      groups.push([a, b]);
      for (const c of cubic) if (c.id > b.id && related(a, c) && related(b, c)) groups.push([a, b, c]);
    }
  });

  const explained = (g: readonly IndexEntry[]) => union(g.map(of));
  const ownCount = (m: IndexEntry, g: readonly IndexEntry[]) => {
    const others = union(g.filter((o) => o !== m).map(of));
    return [...of(m)].filter((k) => !others.has(k)).length;
  };
  const scoreSum = (g: readonly IndexEntry[]) => g.reduce((s, m) => s + (m.score ?? 0), 0);
  const ranked = groups
    .filter((g) => g.every((m) => ownCount(m, g) >= MIN_OWN))
    .map((g) => ({ g, n: explained(g).size, s: scoreSum(g) }))
    .sort((x, y) => y.n - x.n || y.s - x.s || x.g.map((m) => m.id).join().localeCompare(y.g.map((m) => m.id).join()));

  const kept: IndexEntry[][] = [];
  const used = new Set<number>();
  const claimed = new Set<string>();
  for (const { g } of ranked) {
    if (g.some((m) => used.has(m.id))) continue;
    const mine = explained(g);
    if ([...mine].filter((k) => !claimed.has(k)).length < MIN_OWN) continue;
    kept.push(g);
    for (const m of g) used.add(m.id);
    for (const k of mine) claimed.add(k);
  }
  return kept;
}

export function callRows(candidates: readonly IndexEntry[]): CallRow[] {
  const groupOf = new Map<number, IndexEntry[]>();
  for (const g of bonnetGroups(candidates)) for (const m of g) groupOf.set(m.id, g);
  const rows: CallRow[] = [];
  const placed = new Set<IndexEntry[]>();
  for (const c of candidates) {
    const g = groupOf.get(c.id);
    if (!g) {
      rows.push({ kind: "one", candidate: c });
      continue;
    }
    if (placed.has(g)) continue;
    placed.add(g);
    const members = candidates.filter((m) => g.includes(m));
    rows.push({ kind: "group", members, peaks: union(members.map(peakSet)).size });
  }
  return rows;
}

export type KeyedRow = { kind: "one"; candidate: IndexEntry } | { kind: "group"; members: readonly IndexEntry[] };

export function keyedRows(rows: readonly CallRow[]): KeyedRow[] {
  return rows.flatMap((row): KeyedRow[] =>
    row.kind === "one"
      ? [row]
      : [{ kind: "group", members: row.members }, ...row.members.map((candidate) => ({ kind: "one" as const, candidate }))],
  );
}
