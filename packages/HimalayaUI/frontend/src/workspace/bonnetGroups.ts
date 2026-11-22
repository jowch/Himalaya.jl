import type { IndexEntry } from "../api";

import { gaussianCurvature } from "./curvature";

// physics: Bonnet-related D/P/G surfaces share Gaussian curvature κ, so a lattice matches its partner when its κ-equal lattice agrees within K standard errors of the two fits.
const K = 2;
// physics: relative lattice error no fit can beat, set by q binning and detector resolution; it also stands in for fits with too few peaks to estimate their own.
const FLOOR = 0.003;
// physics: bounds how far scattered fits can widen the tolerance; near-aliases such as Ia3d at √(5/2)·a sit inside any tolerance and only MIN_OWN keeps them out.
const CAP = 0.02;
// physics: a κ match alone admits lattice aliases, so every member must explain peaks no other member does; this, not the tolerance, is the alias guard.
const MIN_OWN = 2;

export type CallRow =
  | { kind: "one"; candidate: IndexEntry }
  | { kind: "group"; members: readonly IndexEntry[]; peaks: number };

const peakKey = (q: number) => q.toFixed(6);
const peakSet = (ix: IndexEntry) => new Set(ix.peaks.map((p) => peakKey(p.q_observed)));

// physics: each assigned peak estimates the lattice as a·q_pred/q_obs, so the spread of ln(q_pred/q_obs) over √n is the fit's relative standard error.
function latticeError(ix: IndexEntry): number {
  const logs = ix.peaks.flatMap((p) => {
    const q = ix.predicted_q[p.ratio_position - 1];
    return q !== undefined && q > 0 && p.q_observed > 0 ? [Math.log(q / p.q_observed)] : [];
  });
  const n = logs.length;
  if (n < 3) return FLOOR;
  const mean = logs.reduce((s, x) => s + x, 0) / n;
  const variance = logs.reduce((s, x) => s + (x - mean) ** 2, 0) / (n - 1);
  return Math.hypot(Math.sqrt(variance / n), FLOOR);
}

function related(a: IndexEntry, b: IndexEntry): boolean {
  const ka = gaussianCurvature(a.phase, a.lattice_d);
  const kb = gaussianCurvature(b.phase, b.lattice_d);
  if (ka === null || kb === null || b.lattice_d === null || a.phase === b.phase) return false;
  const tolerance = Math.min(CAP, K * Math.hypot(latticeError(a), latticeError(b)));
  return Math.abs(Math.log(kb / ka)) / 2 <= tolerance;
}

function union(sets: readonly Set<string>[]): Set<string> {
  const out = new Set<string>();
  for (const s of sets) for (const k of s) out.add(k);
  return out;
}

export function bonnetGroups(candidates: readonly IndexEntry[]): IndexEntry[][] {
  const cubic = candidates
    .filter((c) => gaussianCurvature(c.phase, c.lattice_d) !== null)
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
    const fresh = (m: IndexEntry) => {
      const others = union(g.filter((o) => o !== m).map(of));
      return [...of(m)].filter((k) => !others.has(k) && !claimed.has(k)).length;
    };
    if (!g.every((m) => fresh(m) >= MIN_OWN)) continue;
    const mine = explained(g);
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
