/**
 * peak_removed mutator (M2.2). Optimistically removes the peak from the cache
 * and restores on rollback. The DELETE route returns 200 with
 * `{event_id, view_row_id, analysis_inputs_hash}`; onSuccess writes the new
 * hash onto the exposure cache so the StaleIndicesBanner doesn't flash
 * between the optimistic delete and the SSE `post_state` arrival.
 */
import * as api from "../../../api";
import type { Peak, PeakRemoveResponse, Exposure, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export type PeakRemoveInput = { peakId: number };
type PeakRemoveScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const peakRemoveMutator: Mutator<PeakRemoveInput, PeakRemoveScope, PeakRemoveResponse> = {
  kind: "peak_removed",
  onMutate: (p, qc): RollbackContext => {
    // Invariant: optimistic placeholder ids (negative) are an internal cache
    // sentinel and must never be sent over the wire. The TraceViewer click
    // handler skips negative-id peaks, so this should be unreachable — fail
    // loud if anything else routes a negative id here. Without the guard,
    // DELETE /api/peaks/-N would 404, rollback would re-insert the placeholder,
    // and the still-pending add would then resolve to a real peak the user
    // believed they removed (the "ghost peak" desync).
    if (p.peakId < 0) {
      throw new Error(
        `peakRemove invariant: cannot remove optimistic placeholder peak ` +
        `(id=${p.peakId}). The TraceViewer click handler is expected to ` +
        `skip negative-id peaks until the add confirms.`,
      );
    }
    const peaksKey = queryKeys.peaks(p.exposureId);
    const prev = qc.getQueryData<Peak[]>(peaksKey);
    if (prev) {
      // Only manual peaks are removable. The backend's auto_peaks and
      // peak_curations tables are independent SQLite sequences and may share
      // an id, so filtering by id alone could drop the wrong peak. Match on
      // (id, source='manual') to scope to the correct row.
      qc.setQueryData<Peak[]>(peaksKey, prev.filter((pk) =>
        !(pk.id === p.peakId && pk.source === "manual")));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(peaksKey, prev);
      },
    };
  },
  request: (p) => api.removePeak(p.peakId, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // The peak is already removed optimistically. Write the new hash onto the
    // exposure cache so the StaleIndicesBanner doesn't flash before the SSE
    // post_state arrives.
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash: response.analysis_inputs_hash } : old);
  },
  affectsExposurePeaks: () => true,
  // 404 = "already gone" → desired end state. Without this, a 5xx-then-retry
  // can land a 404 on a peak the first attempt already deleted, and rollback
  // re-inserts a phantom row visible until the next refetch.
  // Side note: TanStack does not fire `onSuccess` on a rejected mutationFn,
  // so the `analysis_inputs_hash` write above is skipped on the 404 path.
  // SSE `post_state` will deliver the fresh hash; `useExposureHasPendingPeakOps`
  // masks `StaleIndicesBanner` until then.
  treats404AsSuccess: true,
};
