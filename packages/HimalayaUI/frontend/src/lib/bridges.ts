import type { SampleBridge } from "../api";

function groupOf(list: readonly SampleBridge[], sampleId: number): number | undefined {
  return list.find((r) => r.sample_id === sampleId)?.group_id;
}

function dropSingletons(list: readonly SampleBridge[]): SampleBridge[] {
  const sizes = new Map<number, number>();
  for (const r of list) sizes.set(r.group_id, (sizes.get(r.group_id) ?? 0) + 1);
  return list.filter((r) => (sizes.get(r.group_id) ?? 0) >= 2);
}

export function bridgedWith(list: readonly SampleBridge[], sampleId: number): number[] {
  const g = groupOf(list, sampleId);
  if (g === undefined) return [];
  return list.filter((r) => r.group_id === g && r.sample_id !== sampleId).map((r) => r.sample_id);
}

export function bridgeCounts(list: readonly SampleBridge[]): Map<number, number> {
  const sizes = new Map<number, number>();
  for (const r of list) sizes.set(r.group_id, (sizes.get(r.group_id) ?? 0) + 1);
  return new Map(list.map((r) => [r.sample_id, (sizes.get(r.group_id) ?? 1) - 1]));
}

export function applyBridge(list: readonly SampleBridge[], a: number, b: number): SampleBridge[] {
  if (a === b) return [...list];
  const ga = groupOf(list, a);
  const gb = groupOf(list, b);
  if (ga !== undefined && ga === gb) return [...list];
  if (ga === undefined && gb === undefined) {
    const fresh = Math.min(0, ...list.map((r) => r.group_id)) - 1;
    return [...list, { sample_id: a, group_id: fresh }, { sample_id: b, group_id: fresh }];
  }
  if (ga === undefined) return [...list, { sample_id: a, group_id: gb! }];
  if (gb === undefined) return [...list, { sample_id: b, group_id: ga }];
  return list.map((r) => (r.group_id === gb ? { ...r, group_id: ga } : r));
}

export function applyUnbridge(list: readonly SampleBridge[], sampleId: number): SampleBridge[] {
  return dropSingletons(list.filter((r) => r.sample_id !== sampleId));
}

interface MateCandidate {
  name: string;
  setup_id?: number | null | undefined;
}

function nameKey(name: string): string {
  return name.toLowerCase().replace(/[^a-z0-9]/g, "");
}

function sharedPrefix(a: string, b: string): number {
  let n = 0;
  while (n < a.length && n < b.length && a[n] === b[n]) n++;
  return n;
}

function mateRank(here: MateCandidate, other: MateCandidate): [number, number, number] {
  const sameSetup = here.setup_id != null && other.setup_id === here.setup_id ? 1 : 0;
  const a = nameKey(here.name);
  const b = nameKey(other.name);
  const tier = a === b ? 0 : a.length > 0 && b.length > 0 && (a.includes(b) || b.includes(a)) ? 1 : 2;
  return [sameSetup, tier, -sharedPrefix(a, b)];
}

export function rankBridgeMates<T extends MateCandidate>(here: MateCandidate, others: readonly T[]): T[] {
  const ranked = others.map((s, i) => ({ s, i, r: mateRank(here, s) }));
  ranked.sort((x, y) => x.r[0] - y.r[0] || x.r[1] - y.r[1] || x.r[2] - y.r[2] || x.i - y.i);
  return ranked.map((x) => x.s);
}
