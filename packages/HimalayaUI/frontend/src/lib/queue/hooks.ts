import { useMutationState } from "@tanstack/react-query";
import type { OpKind, OpPayload } from "./types";

/**
 * Op kinds that change the effective peak set for an exposure. Used by
 * `useExposureHasPendingPeakOps` to gate UI elements that should hide or
 * suspend during a peak-curation chain (the stale-indices alert, speculative-snap).
 *
 * `reanalyze_exposure` is included because a manual reanalyze races with
 * curation ops and can produce intermediate states the user shouldn't see.
 */
const PEAK_AFFECTING_KINDS: ReadonlySet<OpKind> = new Set<OpKind>([
  "peak_added",
  "peak_excluded",
  "peak_unexcluded",
  "peak_removed",
  "reanalyze_exposure",
]);

/**
 * True when any peak-affecting mutation is pending for `exposureId`.
 * Returns false when `exposureId` is undefined (no exposure selected).
 */
export function useExposureHasPendingPeakOps(
  exposureId: number | undefined,
): boolean {
  const pending = useMutationState({
    filters: {
      predicate: (m) => {
        if (m.state.status !== "pending") return false;
        if (exposureId === undefined) return false;
        const op = m.state.variables as OpPayload | undefined;
        if (!op) return false;
        if (!PEAK_AFFECTING_KINDS.has(op.kind)) return false;
        return op.exposureId === exposureId;
      },
    },
  });
  return pending.length > 0;
}

