/**
 * peak_added mutator (M2.2). Optimistic effect inserts a manual Peak with a
 * negative placeholder id; onSuccess swaps in the server's authoritative row
 * and updates `exposure.analysis_inputs_hash` so the stale-indices alert does
 * not flash. Indices are not invalidated here — they arrive via SSE
 * `post_state` through `applyRemoteToCache` (replay-without-refetch).
 */
import * as api from "../../../api";
import type { Peak, PeakAddResponse, Exposure, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { nextOptimisticId } from "../optimisticId";
import { peakQTol } from "../peakQTol";
import { replacePlaceholder } from "../replacePlaceholder";
import { stripQueueMetadata } from "../queueMeta";
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
    const { meta, payload: serverPeak } = stripQueueMetadata(response);
    const peaksKey = queryKeys.peaks(p.exposureId);
    qc.setQueryData<Peak[]>(peaksKey, (old) =>
      // Dedup is scoped to MANUAL peaks: auto_peaks.id and peak_curations.id
      // share a wire namespace, so an auto peak with the same id must survive.
      replacePlaceholder(
        old ?? [],
        serverPeak,
        (pk) => pk.source !== "auto"
             && Math.abs(pk.q - p.q) < peakQTol(p.q)
             && pk.exposure_id === p.exposureId,
        { isDuplicate: (pk) => pk.source === "manual" && pk.id === serverPeak.id },
      ),
    );
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash: meta.analysis_inputs_hash ?? null } : old);
  },
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    const peakId = payload.peak_curation_id as number | undefined;
    if (peakId === undefined) return undefined;
    return {
      ...base,
      id: peakId,
      exposure_id: remote.entity_id,
      q: payload.q as number,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
      view_row_id: peakId,
    } as PeakAddResponse;
  },
  affectsExposurePeaks: () => true,
};
