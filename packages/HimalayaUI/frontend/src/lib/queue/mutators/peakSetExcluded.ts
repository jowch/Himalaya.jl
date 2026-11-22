import * as api from "../../../api";
import type { Peak, PeakUpdatedResponse, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpKind, RollbackContext } from "../types";
import { refreshPeaks } from "../refreshPeaks";
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
    request: (p) => api.setAutoPeakExcluded(p.peakId, excluded, buildAuthOpts(p)),
    onSuccess: (p, response, qc) => {
      const peaksKey = queryKeys.peaks(p.exposureId);
      const { payload: peakFields } = stripQueueMetadata(response);
      qc.setQueryData<Peak[]>(peaksKey, (old) =>
        (old ?? []).map((pk) =>
          pk.id === peakFields.id && pk.source === "auto"
            ? ({ ...pk, ...peakFields } as Peak) : pk,
        ),
      );
      refreshPeaks(qc, p.exposureId);
    },
    synthesizeFromSse: (remote, base) => {
      const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
      if (payload.auto_peak_id === undefined) return undefined;
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
