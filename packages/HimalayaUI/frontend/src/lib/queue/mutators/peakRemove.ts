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

export const peakRemoveMutator: Mutator<OpPayload<PeakRemoveInput>, PeakRemoveResponse> = {
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
  onSuccess: (raw, response, qc) => {
    const p = flat(raw);
    // The peak is already removed optimistically. Write the new hash onto the
    // exposure cache so the StaleIndicesBanner doesn't flash before the SSE
    // post_state arrives.
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash: response.analysis_inputs_hash } : old);
  },
  affectsExposurePeaks: () => true,
};
