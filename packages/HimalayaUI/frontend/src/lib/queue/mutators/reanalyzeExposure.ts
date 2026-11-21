/**
 * reanalyze_exposure mutator (M2.5). No UI sends it today; the registry keeps
 * it so a persisted op still rehydrates. There is no real optimistic effect — the server runs
 * analysis and the resulting indices arrive via SSE post_state on the
 * analyze_run frame (handled by applyRemoteToCache.ts). The onMutate returns
 * a no-op restore.
 *
 * Migrated from a legacy useMutation so it composes with M0 fast-skip
 * (no double round-trip when peak ops just produced the latest hash) and
 * preserves FIFO ordering with peak ops via MutationCache iteration order.
 */
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
  // Null optimistic effect: no cache write, restore is a no-op.
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.reanalyzeExposure(p.exposureId, buildAuthOpts(p)),
  onSuccess: () => {},
};
