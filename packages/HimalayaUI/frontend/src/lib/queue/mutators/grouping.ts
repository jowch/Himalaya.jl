import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { Load, LoadSample } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";          // src/lib/authOpts.ts — the shared helper
import type { Mutator, RollbackContext } from "../types";

// The shared scope every grouping mutator carries (experiment + identity).
// merge/split overlap input and scope; see Task 5.
interface BaseScope {
  username: string | undefined;
  clientId: string;
  experimentId: number;
}

// authOpts(username, clientId, clientOpId) — same idiom as trivial.ts:29-35.
function buildAuthOpts(p: {
  username?: string | undefined;
  clientId?: string | undefined;
  clientOpId?: string | undefined;
}): api.AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/** Invalidate the loads roll-up for an experiment. Called from EVERY grouping
 *  mutator's onSuccess: the replay coordinator's own-op path resolves the
 *  deferred + aborts the HTTP request and NEVER calls applyRemoteToCache
 *  (replayCoordinator.ts case-1), so the SSE arm alone would not reconcile an
 *  own-op edit — the saveSeries.ts `series_save` precedent invalidates in
 *  onSuccess for exactly this reason. E1's applyRemoteToCache structural arms
 *  (E1-owned; E2 verifies in Task 7) cover FOREIGN tabs. */
function invalidateLoads(qc: QueryClient, experimentId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.loads(experimentId) });
  qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) });
}

/** Snapshot the loads cache, run `mutate`, return a restore closure. Shared by
 *  every grouping mutator so rollback is uniform. */
function patchLoads(
  qc: QueryClient, experimentId: number, mutate: (loads: Load[]) => Load[],
): RollbackContext {
  const key = queryKeys.loads(experimentId);
  const prev = qc.getQueryData<Load[]>(key);
  if (prev) qc.setQueryData<Load[]>(key, mutate(structuredClone(prev)));
  return { restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); } };
}

// ---------------------------------------------------------------------------
// move_exposure  (single-entity → exposure_moved)
// Entity ids in Input, NOT scope — one hook instance per experiment, not per row.
// ---------------------------------------------------------------------------
type MoveExposureInput = { exposureId: number; sampleId: number };
type MoveExposureScope = BaseScope;

export const moveExposureMutator: Mutator<MoveExposureInput, MoveExposureScope, api.Exposure> = {
  kind: "move_exposure",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      let moved: LoadSample["exposures"][number] | undefined;
      for (const ld of loads) {
        for (const s of ld.samples) {
          const i = s.exposures.findIndex((e) => e.id === p.exposureId); // §8.8: exposure key is `.id`
          if (i >= 0) { moved = s.exposures.splice(i, 1)[0]; break; }
        }
        if (moved) break;
      }
      if (moved) {
        for (const ld of loads) {
          const dest = ld.samples.find((s) => s.sample_id === p.sampleId);
          if (dest) { dest.exposures.push(moved); break; }
        }
      }
      return loads;
    }),
  request: (p) => api.moveExposure(p.exposureId, p.sampleId, buildAuthOpts(p)),
  // Own-op reconcile: the surgical onMutate reflects the end state, but a move
  // can flip a flag (re-derived server-side), and replay never calls
  // applyRemoteToCache for own ops — so invalidate loads(id) to refetch.
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// rename_sample  (single-entity → sample_renamed)
// Entity id in Input, NOT scope — one hook instance per experiment, not per row.
// ---------------------------------------------------------------------------
type RenameSampleInput = { sampleId: number; name: string };
type RenameSampleScope = BaseScope;

export const renameSampleMutator: Mutator<RenameSampleInput, RenameSampleScope, api.Sample> = {
  kind: "rename_sample",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s) { s.name = p.name; s.name_source = "user"; break; }
      }
      return loads;
    }),
  request: (p) => api.renameSample(p.sampleId, p.name, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // `response.name` may be a synthesized PARTIAL on the own-op SSE-wins path
    // (replayCoordinator resolves the deferred with a frame-derived stub, not
    // the HTTP body — feedback_queue_sse_wins_partial), so guard it before
    // splicing; then invalidate to re-derive flags / refresh the flat samples.
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s && response?.name) s.name = response.name;
      }
      return loads;
    });
    invalidateLoads(qc, p.experimentId);
  },
};
