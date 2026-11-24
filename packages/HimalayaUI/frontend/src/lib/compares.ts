import type { CompareChange, CompareState, SampleBackground } from "../api";

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
  const others = state.subtractions.filter((t) => t.sample_id !== sampleId);
  const background = state.backgrounds?.find((b) => b.sample_id === sampleId);
  const bgOthers = (state.backgrounds ?? []).filter((b) => b.sample_id !== sampleId);
  if (background) {
    const wanted = change.subtraction;
    const nextBackground: SampleBackground =
      wanted === undefined ? background
      : wanted === "default" ? { ...background, ref_sample_id: background.default_ref_sample_id, source: "default" }
      : { ...background, ref_sample_id: wanted, source: "user" };
    return {
      compares: [...kept, ...added],
      subtractions: typeof wanted === "number" ? [...others, { sample_id: sampleId, ref_sample_id: wanted }] : wanted === undefined ? state.subtractions : others,
      backgrounds: [...bgOthers, nextBackground],
    };
  }
  const current = state.subtractions.find((t) => t.sample_id === sampleId)?.ref_sample_id;
  const stale = current !== undefined && !before.has(current) && added.length + removing.length > 0;
  const requested = change.subtraction === "default" ? null : change.subtraction;
  const nextSub = "subtraction" in change ? requested : current !== undefined && (removing.includes(current) || stale) ? null : current;
  return {
    compares: [...kept, ...added],
    subtractions: nextSub === null || nextSub === undefined ? others : [...others, { sample_id: sampleId, ref_sample_id: nextSub }],
    ...(state.backgrounds ? { backgrounds: state.backgrounds } : {}),
  };
}

export function comparesOf(state: CompareState, sampleId: number): { refs: number[]; addedBy: ReadonlyMap<number, string | null>; subtraction: number | undefined; background: SampleBackground | undefined; chosen: number | null | undefined } {
  const mine = state.compares.filter((c) => c.sample_id === sampleId);
  const background = state.backgrounds?.find((b) => b.sample_id === sampleId);
  const subtraction = state.subtractions.find((t) => t.sample_id === sampleId)?.ref_sample_id;
  return {
    refs: mine.map((c) => c.ref_sample_id),
    addedBy: new Map(mine.map((c) => [c.ref_sample_id, c.added_by])),
    subtraction,
    background,
    chosen: background ? (background.source === "default" ? undefined : background.ref_sample_id) : subtraction,
  };
}
