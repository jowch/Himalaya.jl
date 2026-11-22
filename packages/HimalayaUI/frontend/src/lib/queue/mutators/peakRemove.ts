import * as api from "../../../api";
import type { Peak, PeakRemoveResponse, AuthOpts } from "../../../api";
import { cancelPeakRefetch, refreshPeaks } from "../refreshPeaks";
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
    if (p.peakId < 0) {
      throw new Error(
        `peakRemove invariant: cannot remove optimistic placeholder peak ` +
        `(id=${p.peakId}). The TraceViewer click handler is expected to ` +
        `skip negative-id peaks until the add confirms.`,
      );
    }
    cancelPeakRefetch(qc, p.exposureId);
    const peaksKey = queryKeys.peaks(p.exposureId);
    const prev = qc.getQueryData<Peak[]>(peaksKey);
    if (prev) {
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
  onSuccess: (p, _response, qc) => refreshPeaks(qc, p.exposureId),
  treats404AsSuccess: true,
};
