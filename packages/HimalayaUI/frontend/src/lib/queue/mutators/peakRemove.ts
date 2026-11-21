/**
 * peak_removed mutator (M2.2). Optimistically removes the peak from the cache
 * and restores on rollback. Indices arrive via the SSE `post_state`.
 */
import * as api from "../../../api";
import type { Peak, PeakRemoveResponse, AuthOpts } from "../../../api";
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
    // handler skips negative-id peaks, so this should be unreachable. If a
    // future regression at the click layer routes a negative id here, the
    // throw surfaces as a console error (the thrown Error has no `status`,
    // so it doesn't route through the validation toast or retry path; the
    // mutation simply settles in error state and the user sees no UI change).
    // That's intentional: the user did nothing wrong, and the desync this
    // guard prevents — DELETE /api/peaks/manual/-N → 404 → rollback re-inserts the
    // placeholder while the still-pending add resolves into a real "ghost
    // peak" — should fail at the dev/test layer, not surface as a toast.
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
  request: (p) => api.removeManualPeak(p.peakId, buildAuthOpts(p)),
  onSuccess: () => {},
  // 404 = "already gone" → desired end state. Without this, a 5xx-then-retry
  // can land a 404 on a peak the first attempt already deleted, and rollback
  // re-inserts a phantom row visible until the next refetch.
  treats404AsSuccess: true,
};
