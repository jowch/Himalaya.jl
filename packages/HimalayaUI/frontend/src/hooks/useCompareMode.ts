/**
 * `CompareMode` — the tagged union describing what the Compare surface is
 * currently doing (spec §3 / Compare UX refinement).
 *
 * The `useCompareMode` hook derives a `CompareMode` from `activeDraft` + the
 * loaded `Comparison` (Compare UX C-5). The leaf components (SavePill, etc.)
 * only need the type, so it was split out here to keep #139 self-contained
 * and free of `ActiveDraft`-shape coupling.
 */
import { useMemo } from "react";
import { useAppState } from "../state";
import type { Comparison } from "../api";

export type CompareMode =
  | { kind: "viewing" }
  | { kind: "viewing-stale"; previousHash: string }
  | { kind: "editing-mine"; draftId: number }
  | { kind: "editing-as-fork-of"; parentId: number }
  | { kind: "creating-blank" }
  | { kind: "creating-from-fork"; parentId: number };

export function useCompareMode(opts: {
  comparison: Comparison | undefined;
  currentUserId: number | undefined;
  /** Set when a foreign SSE event drifted content_hash since last read. */
  staleAgainstHash?: string | undefined;
}): CompareMode {
  const draft = useAppState((s) => s.activeDraft);

  return useMemo<CompareMode>(() => {
    if (draft === null) {
      if (opts.staleAgainstHash !== undefined) {
        return { kind: "viewing-stale", previousHash: opts.staleAgainstHash };
      }
      return { kind: "viewing" };
    }
    if (draft.id === undefined) {
      if (draft.forkedFromId !== undefined) {
        return { kind: "creating-from-fork", parentId: draft.forkedFromId };
      }
      return { kind: "creating-blank" };
    }
    const author = opts.comparison?.created_by ?? null;
    const isAuthor =
      author !== null &&
      opts.currentUserId !== undefined &&
      opts.currentUserId === author;
    if (isAuthor) return { kind: "editing-mine", draftId: draft.id };
    return { kind: "editing-as-fork-of", parentId: draft.id };
  }, [draft, opts.currentUserId, opts.comparison?.created_by, opts.staleAgainstHash]);
}
