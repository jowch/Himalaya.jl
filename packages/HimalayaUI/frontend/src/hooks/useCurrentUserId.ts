/**
 * useCurrentUserId — resolve the current user's numeric id from the
 * Zustand `username` against the cached `listUsers` query.
 *
 * Returns `undefined` until resolution lands (or if no user is selected).
 * Mirrors the inline copies that lived in `NeedsReviewBadge.tsx` and
 * `ComparisonPicker.tsx`; promoted to a shared hook in Phase 11 because
 * the Edit-vs-Fork affordance also needs it.
 */
import { useQuery } from "@tanstack/react-query";
import { useAppState } from "../state";
import * as api from "../api";

export function useCurrentUserId(): number | undefined {
  const username = useAppState((s) => s.username);
  const usersQ = useQuery({
    queryKey: ["users"] as const,
    queryFn: () => api.listUsers(),
    enabled: username !== undefined,
  });
  if (username === undefined) return undefined;
  const u = (usersQ.data ?? []).find((x) => x.username === username);
  return u?.id;
}
