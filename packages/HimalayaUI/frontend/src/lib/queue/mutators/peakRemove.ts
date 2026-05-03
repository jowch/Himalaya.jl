/**
 * peak_removed mutator (M2.2). Optimistically removes the peak from the cache
 * and restores on rollback. The DELETE route returns 204; the post-state
 * (exposure hash + indices) arrives via SSE `post_state` through
 * `applyRemoteToCache` — no work needed in onSuccess.
 */
import * as api from "../../../api";
import type { Peak, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpPayload, RollbackContext } from "../types";

export type PeakRemoveInput = { peakId: number };
type PeakRemoveScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};
type Flat = OpPayload<PeakRemoveInput> & PeakRemoveScope & PeakRemoveInput;
const flat = (p: OpPayload<PeakRemoveInput>): Flat => p as unknown as Flat;

function buildAuthOpts(p: Flat): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const peakRemoveMutator: Mutator<OpPayload<PeakRemoveInput>, void> = {
  kind: "peak_removed",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat(raw);
    const peaksKey = queryKeys.peaks(p.exposureId);
    const prev = qc.getQueryData<Peak[]>(peaksKey);
    if (prev) {
      qc.setQueryData<Peak[]>(peaksKey, prev.filter((pk) => pk.id !== p.peakId));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(peaksKey, prev);
      },
    };
  },
  request: (raw) => {
    const p = flat(raw);
    return api.removePeak(p.peakId, buildAuthOpts(p));
  },
  onSuccess: () => {
    // No-op: the peak is already removed optimistically. SSE post_state will
    // refresh the exposure hash + indices via applyRemoteToCache.
  },
  affectsExposurePeaks: () => true,
};
