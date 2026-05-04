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
      const {
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        event_id, view_row_id, analysis_inputs_hash, client_op_id,
        ...peakFields
      } = response as PeakUpdatedResponse & { client_op_id?: string };
      qc.setQueryData<Peak[]>(peaksKey, (old) =>
        (old ?? []).map((pk) =>
          pk.id === peakFields.id && pk.source === "auto"
            ? ({ ...pk, ...peakFields } as Peak) : pk,
        ),
      );
      qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
        old ? { ...old, analysis_inputs_hash: response.analysis_inputs_hash } : old);
    },
    affectsExposurePeaks: () => true,
  };
}

export const peakExcludeMutator = makeMutator("peak_excluded", true);
export const peakUnexcludeMutator = makeMutator("peak_unexcluded", false);
