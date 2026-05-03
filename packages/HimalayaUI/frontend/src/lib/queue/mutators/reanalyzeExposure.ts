/**
 * reanalyze_exposure mutator (M2.5). Manual "Re-analyze" button on
 * StaleIndicesBanner. There is no real optimistic effect — the server runs
 * analysis and the resulting indices/groups arrive via SSE post_state on the
 * analyze_run frame (or are written to the cache here as defense-in-depth
 * when SSE is unavailable). The onMutate returns a no-op restore.
 *
 * Migrated from a legacy useMutation so it composes with M0 fast-skip
 * (no double round-trip when peak ops just produced the latest hash) and
 * preserves FIFO ordering with peak ops via MutationCache iteration order.
 */
import * as api from "../../../api";
import type { AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpPayload, RollbackContext } from "../types";

export type ReanalyzeExposureInput = Record<string, never>;

type ReanalyzeExposureScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

type Flat = OpPayload<ReanalyzeExposureInput> & ReanalyzeExposureScope;
const flat = (p: OpPayload<ReanalyzeExposureInput>): Flat =>
  p as unknown as Flat;

function buildAuthOpts(p: Flat): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

type ReanalyzeResponse = { id: number; analyzed: boolean };

export const reanalyzeExposureMutator: Mutator<
  OpPayload<ReanalyzeExposureInput>,
  ReanalyzeResponse
> = {
  kind: "reanalyze_exposure",
  // Null optimistic effect: no cache write, restore is a no-op.
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (raw) => {
    const p = flat(raw);
    return api.reanalyzeExposure(p.exposureId, buildAuthOpts(p));
  },
  // The response is purely advisory ({ id, analyzed }); cross-tab indices
  // arrive via SSE post_state on the analyze_run frame. Invalidate scoped
  // queries to ensure a fresh fetch lands even when SSE is offline.
  onSuccess: (raw, _response, qc) => {
    const p = flat(raw);
    qc.invalidateQueries({ queryKey: queryKeys.peaks(p.exposureId) });
    qc.invalidateQueries({ queryKey: queryKeys.indices(p.exposureId) });
    qc.invalidateQueries({ queryKey: queryKeys.groups(p.exposureId) });
  },
  // Marks this op as peak-affecting so StaleIndicesBanner / speculative-snap
  // gating treats an in-flight reanalyze the same way as a peak curation.
  affectsExposurePeaks: () => true,
};
