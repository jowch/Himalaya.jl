import * as api from "../../../api";
import type { AuthOpts, ReanalyzeResponse } from "../../../api";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export type ReanalyzeExposureInput = Record<string, never>;

type ReanalyzeExposureScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const reanalyzeExposureMutator: Mutator<
  ReanalyzeExposureInput,
  ReanalyzeExposureScope,
  ReanalyzeResponse
> = {
  kind: "reanalyze_exposure",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.reanalyzeExposure(p.exposureId, buildAuthOpts(p)),
  onSuccess: () => {},
};
