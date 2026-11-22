import type { CompareChange, CompareState } from "../api";

export const EMPTY_COMPARES: CompareState = { compares: [], subtractions: [] };

export function applyCompareChange(state: CompareState, sampleId: number, change: CompareChange, username: string | null): CompareState {
  const remove = new Set(change.remove ?? []);
  const before = new Set(state.compares.filter((c) => c.sample_id === sampleId).map((c) => c.ref_sample_id));
  const removing = [...remove].filter((r) => before.has(r));
  const kept = state.compares.filter((c) => !(c.sample_id === sampleId && remove.has(c.ref_sample_id)));
  const have = new Set(kept.filter((c) => c.sample_id === sampleId).map((c) => c.ref_sample_id));
  const added = [...new Set(change.add ?? [])]
    .filter((r) => !have.has(r) && !remove.has(r))
    .map((ref_sample_id) => ({ sample_id: sampleId, ref_sample_id, added_by: username }));
  const current = state.subtractions.find((t) => t.sample_id === sampleId)?.ref_sample_id;
  const stale = current !== undefined && !before.has(current) && added.length + removing.length > 0;
  const nextSub = "subtraction" in change ? change.subtraction : current !== undefined && (removing.includes(current) || stale) ? null : current;
  const others = state.subtractions.filter((t) => t.sample_id !== sampleId);
  return {
    compares: [...kept, ...added],
    subtractions: nextSub === null || nextSub === undefined ? others : [...others, { sample_id: sampleId, ref_sample_id: nextSub }],
  };
}

export function comparesOf(state: CompareState, sampleId: number): { refs: number[]; addedBy: ReadonlyMap<number, string | null>; subtraction: number | undefined } {
  const mine = state.compares.filter((c) => c.sample_id === sampleId);
  return {
    refs: mine.map((c) => c.ref_sample_id),
    addedBy: new Map(mine.map((c) => [c.ref_sample_id, c.added_by])),
    subtraction: state.subtractions.find((t) => t.sample_id === sampleId)?.ref_sample_id,
  };
}
