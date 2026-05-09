// prefetchMembers.ts — warm the four cache keys computeMemberSnapshot reads.
//
// Two call sites (handleSave in ComparePageEdit, buildOverwritePayload in
// ConflictModal) both need the same four-key prefetch before computing
// snapshots — without it, never-visited members land with
// analysis_inputs_hash = "" and the server marks them stale on the next
// view fold (issue #49). Extracted here so the call sites can't drift
// (#73-class drift seam, tracked as #93).
//
// Per-key cold detection: each of the four keys is checked independently,
// so an exposure with three warm keys and one cold key only refetches the
// missing one. Cheaper than the previous all-or-nothing per-exposure check
// AND more honest about what's actually being warmed.
import type { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../queries";
import { getExposure, listPeaks, listIndices, listGroups } from "../../api";

export async function prefetchColdMembers(
  exposureIds: ReadonlyArray<number>,
  qc: QueryClient,
): Promise<void> {
  const fetches: Promise<unknown>[] = [];
  for (const id of exposureIds) {
    if (qc.getQueryData(queryKeys.exposure(id)) === undefined) {
      fetches.push(qc.fetchQuery({
        queryKey: queryKeys.exposure(id),
        queryFn: () => getExposure(id),
      }));
    }
    if (qc.getQueryData(queryKeys.peaks(id)) === undefined) {
      fetches.push(qc.fetchQuery({
        queryKey: queryKeys.peaks(id),
        queryFn: () => listPeaks(id),
      }));
    }
    if (qc.getQueryData(queryKeys.indices(id)) === undefined) {
      fetches.push(qc.fetchQuery({
        queryKey: queryKeys.indices(id),
        queryFn: () => listIndices(id),
      }));
    }
    if (qc.getQueryData(queryKeys.groups(id)) === undefined) {
      fetches.push(qc.fetchQuery({
        queryKey: queryKeys.groups(id),
        queryFn: () => listGroups(id),
      }));
    }
  }
  if (fetches.length === 0) return;
  await Promise.all(fetches);
}
