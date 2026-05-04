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
import { nextOptimisticId } from "../optimisticId";
import { peakQTol } from "../peakQTol";
import type { Mutator, RollbackContext } from "../types";

export type PeakAddInput = { q: number };
type PeakAddScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const peakAddMutator: Mutator<PeakAddInput, PeakAddScope, PeakAddResponse> = {
  kind: "peak_added",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.addPeak(p.exposureId, p.q, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // Strip queue-framework metadata before treating the response as a Peak.
    // The route returns a flat `Peak & {event_id, view_row_id, analysis_inputs_hash}`;
    // the cache holds Peak[], so we don't want event_id/etc. polluting it.
    const { event_id: _e, view_row_id: _v, analysis_inputs_hash, ...serverPeak } = response;
    void _e; void _v;
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
            && Math.abs(pk.q - p.q) < peakQTol(p.q)
            && pk.exposure_id === p.exposureId) {
          if (!seen.has(serverPeak.id)) {
            next.push(serverPeak);
            seen.add(serverPeak.id);
          }
          replaced = true;
          continue;
        }
        if (seen.has(pk.id)) continue;
        next.push(pk);
        seen.add(pk.id);
      }
      if (!replaced && !seen.has(serverPeak.id)) next.push(serverPeak);
      return next;
    });
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash } : old);
  },
  affectsExposurePeaks: () => true,
};
