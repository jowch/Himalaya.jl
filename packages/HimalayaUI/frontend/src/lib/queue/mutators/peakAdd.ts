/**
 * peak_added mutator (M2.2). Optimistic effect inserts a manual Peak with a
 * negative placeholder id; onSuccess swaps in the server's authoritative row
 * and updates `exposure.analysis_inputs_hash` so `StaleIndicesBanner` does
 * not flash. Indices are not invalidated here — they arrive via SSE
 * `post_state` through `applyRemoteToCache` (replay-without-refetch).
 */
import * as api from "../../../api";
import type { Peak, PeakAddResponse, Exposure, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpPayload, RollbackContext } from "../types";

export type PeakAddInput = { q: number };
type PeakAddScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};
type Flat = OpPayload<PeakAddInput> & PeakAddScope & PeakAddInput;
const flat = (p: OpPayload<PeakAddInput>): Flat => p as unknown as Flat;

let optimisticCounter = 0;
function nextOptimisticId(): number {
  optimisticCounter += 1;
  return -(Date.now() * 1000 + optimisticCounter);
}

function buildAuthOpts(p: Flat): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const peakAddMutator: Mutator<OpPayload<PeakAddInput>, PeakAddResponse> = {
  kind: "peak_added",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat(raw);
    const peaksKey = queryKeys.peaks(p.exposureId);
    const prev = qc.getQueryData<Peak[]>(peaksKey);
    const placeholder: Peak = {
      id: nextOptimisticId(),
      exposure_id: p.exposureId,
      q: p.q,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
    };
    qc.setQueryData<Peak[]>(peaksKey, [...(prev ?? []), placeholder]);
    return {
      restore: () => {
        if (prev === undefined) qc.removeQueries({ queryKey: peaksKey, exact: true });
        else qc.setQueryData(peaksKey, prev);
      },
    };
  },
  request: (raw) => {
    const p = flat(raw);
    return api.addPeak(p.exposureId, p.q, buildAuthOpts(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flat(raw);
    const peaksKey = queryKeys.peaks(p.exposureId);
    qc.setQueryData<Peak[]>(peaksKey, (old) => {
      const list = old ?? [];
      // Drop the most recent negative-id placeholder for this q (within tol),
      // and dedupe against any concurrent SSE that already inserted the row.
      const seen = new Set<number>();
      const next: Peak[] = [];
      let replaced = false;
      for (const pk of list) {
        if (pk.id < 0 && !replaced
            && Math.abs(pk.q - p.q) < 1e-9
            && pk.exposure_id === p.exposureId) {
          if (!seen.has(response.peak.id)) {
            next.push(response.peak);
            seen.add(response.peak.id);
          }
          replaced = true;
          continue;
        }
        if (seen.has(pk.id)) continue;
        next.push(pk);
        seen.add(pk.id);
      }
      if (!replaced && !seen.has(response.peak.id)) next.push(response.peak);
      return next;
    });
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash: response.analysis_inputs_hash } : old);
  },
  affectsExposurePeaks: () => true,
};
