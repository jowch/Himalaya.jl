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
