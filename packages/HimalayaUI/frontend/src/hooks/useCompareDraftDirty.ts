import { useAppState } from "../state";

/**
 * Returns `true` whenever there is any active draft. The simplification —
 * "presence of a draft = dirty" — is intentional. Drafts only exist as
 * the result of an actual edit gesture (per spec §3, baseHash capture
 * happens on first edit), so this is sufficient.
 */
export function useCompareDraftDirty(): boolean {
  return useAppState((s) => s.activeDraft !== null);
}
