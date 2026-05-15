/**
 * peak_excluded / peak_unexcluded mutators (M2.2). Backend uses two distinct
 * event kinds depending on the boolean target, so we ship two mutators with
 * the same payload + cache logic and route between them at the hook layer
 * based on `excluded`. Optimistically flips the peak's `excluded` field;
 * onSuccess writes the server peak (authoritative) and updates the exposure
 * `analysis_inputs_hash` so `StaleIndicesBanner` does not flash.
 */
import * as api from "../../../api";
import type { Peak, PeakUpdatedResponse, Exposure, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpKind, RollbackContext } from "../types";
import { stripQueueMetadata } from "../queueMeta";

export type PeakSetExcludedInput = { peakId: number };
type PeakSetExcludedScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

function makeMutator(
  kind: OpKind,
  excluded: boolean,
): Mutator<PeakSetExcludedInput, PeakSetExcludedScope, PeakUpdatedResponse> {
  return {
    kind,
    onMutate: (p, qc): RollbackContext => {
      const peaksKey = queryKeys.peaks(p.exposureId);
      const prev = qc.getQueryData<Peak[]>(peaksKey);
      if (prev) {
        // Only auto peaks are excludable. Match (id, source='auto') so a
        // manual peak that happens to share the id (independent SQLite
        // sequences) does not get its excluded flag flipped.
        qc.setQueryData<Peak[]>(peaksKey, prev.map((pk) =>
          pk.id === p.peakId && pk.source === "auto" ? { ...pk, excluded } : pk,
        ));
      }
      return {
        restore: () => {
          if (prev !== undefined) qc.setQueryData(peaksKey, prev);
        },
      };
    },
    request: (p) => api.setPeakExcluded(p.peakId, excluded, buildAuthOpts(p)),
    onSuccess: (p, response, qc) => {
      const peaksKey = queryKeys.peaks(p.exposureId);
      // Strip event metadata + queue plumbing fields. Merge (not replace) so
      // SSE-wins synthetic responses — which omit intensity/prominence/
      // sharpness — preserve those fields from the optimistic row. HTTP-wins
      // responses include the full Peak shape, so merge is equivalent to
      // replace there.
      const { meta, payload: peakFields } = stripQueueMetadata(response);
      qc.setQueryData<Peak[]>(peaksKey, (old) =>
        (old ?? []).map((pk) =>
          pk.id === peakFields.id && pk.source === "auto"
            ? ({ ...pk, ...peakFields } as Peak) : pk,
        ),
      );
      qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
        old ? { ...old, analysis_inputs_hash: meta.analysis_inputs_hash ?? null } : old);
    },
    affectsExposurePeaks: () => true,
    synthesizeFromSse: (remote, base) => {
      const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
      if (payload.auto_peak_id === undefined) return undefined;
      // intensity/prominence/sharpness/exposure_id are intentionally absent —
      // they aren't carried on the SSE wire. The cast goes through `unknown`
      // because PeakUpdatedResponse extends Peak (which requires those fields),
      // and `as PeakUpdatedResponse` alone would be rejected by tsc-strict for
      // insufficient overlap. Safe at runtime: peakSetExcluded.onSuccess uses
      // spread merge `{ ...pk, ...peakFields }` into the existing cached Peak,
      // so the optimistic row's detection values survive when the synth omits
      // them. The legacy synth in replayCoordinator relied on the same merge —
      // preserve the contract.
      return {
        ...base,
        id: payload.auto_peak_id as number,
        q: payload.q as number,
        source: "auto",
        excluded,
      } as unknown as PeakUpdatedResponse;
    },
  };
}

export const peakExcludeMutator = makeMutator("peak_excluded", true);
export const peakUnexcludeMutator = makeMutator("peak_unexcluded", false);
