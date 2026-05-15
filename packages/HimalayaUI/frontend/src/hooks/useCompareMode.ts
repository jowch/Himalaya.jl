/**
 * `CompareMode` — the tagged union describing what the Compare surface is
 * currently doing (spec §3 / Compare UX refinement).
 *
 * Only the *type* lives here for now. The `useCompareMode` hook that derives
 * a `CompareMode` from `activeDraft` + the loaded `Comparison` is Task C-5,
 * landed by the Phase C integration issue (#142) — it appends to this file.
 * The leaf components (SavePill, etc.) only need the type, so it is split out
 * here to keep #139 self-contained and free of `ActiveDraft`-shape coupling.
 */
export type CompareMode =
  | { kind: "viewing" }
  | { kind: "viewing-stale"; previousHash: string }
  | { kind: "editing-mine"; draftId: number }
  | { kind: "editing-as-fork-of"; parentId: number }
  | { kind: "creating-blank" }
  | { kind: "creating-from-fork"; parentId: number };
