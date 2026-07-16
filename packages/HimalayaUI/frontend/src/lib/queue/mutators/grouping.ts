import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { Load, LoadSample } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";          // src/lib/authOpts.ts — the shared helper
import type { Mutator, RollbackContext } from "../types";
import { nextOptimisticId } from "../optimisticId";   // negative placeholder id

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
 *  (E1-owned; E2 verifies in Task 7) cover FOREIGN tabs.
 *
 *  corpusSamples is included to match the foreign arms in applyRemoteToCache
 *  (sample_created / sample_split / grouping_flag_dismissed / sample_renamed all
 *  invalidate corpusSamples) so the contact sheet refreshes on own-tab edits. */
function invalidateLoads(qc: QueryClient, experimentId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.loads(experimentId) });
  qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) });
  qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
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
    const newName = response?.name;
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s && newName) s.name = newName;
      }
      return loads;
    });
    // Mirror the foreign arm's surgical patch (applyRemoteToCache sample_renamed):
    // update the sample-entity cache so the originating tab's sample-detail page
    // converges without waiting for the background refetch.
    if (newName !== undefined) {
      qc.setQueryData<api.Sample>(queryKeys.sample(p.sampleId), (old) =>
        old ? { ...old, name: newName } : old);
    }
    invalidateLoads(qc, p.experimentId);
  },
};

// ---------------------------------------------------------------------------
// merge_samples  (orchestration; NO sample_merged event — the backend fans it
// out as exposure_moved frames. Optimistic tree edit, invalidate-only confirm.)
// ---------------------------------------------------------------------------
type MergeSamplesInput = { loserId: number; survivorId: number };
type MergeSamplesScope = BaseScope;

export const mergeSamplesMutator: Mutator<MergeSamplesInput, MergeSamplesScope, api.MergeSamplesResponse> = {
  kind: "merge_samples",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      // Gather the loser's exposures, drop the loser sample, then append those
      // exposures onto the survivor (preserving order).
      let loserExps: LoadSample["exposures"] = [];
      for (const ld of loads) {
        const li = ld.samples.findIndex((s) => s.sample_id === p.loserId);
        if (li >= 0) { loserExps = ld.samples[li]!.exposures; ld.samples.splice(li, 1); break; }
      }
      for (const ld of loads) {
        const surv = ld.samples.find((s) => s.sample_id === p.survivorId);
        if (surv) {
          surv.exposures.push(...loserExps);
          surv.grouping_source = "user_merged";
          break;
        }
      }
      return loads;
    }),
  request: (p) => api.mergeSamples(p.loserId, p.survivorId, buildAuthOpts(p)),
  // Invalidate-only confirm: the re-derived flags aren't response-derivable, and
  // the own-op path never calls applyRemoteToCache — so refetch loads(id) here.
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// split_sample  (orchestration → sample_created + sample_split; invalidate-only)
// sampleId (source) in Input, NOT scope — one hook per experiment, not per row.
// ---------------------------------------------------------------------------
type SplitSampleInput = { sampleId: number; exposureIds: number[]; name: string };
type SplitSampleScope = BaseScope;

export const splitSampleMutator: Mutator<SplitSampleInput, SplitSampleScope, api.SplitSampleResponse> = {
  kind: "split_sample",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      const ids = new Set(p.exposureIds);
      for (const ld of loads) {
        const src = ld.samples.find((s) => s.sample_id === p.sampleId);
        if (!src) continue;
        const moved = src.exposures.filter((e) => ids.has(e.id));        // §8.8: `.id`
        src.exposures = src.exposures.filter((e) => !ids.has(e.id));
        const created: LoadSample = {
          sample_id: nextOptimisticId(), // NEGATIVE placeholder until SSE confirms
          name: p.name, slot_index: src.slot_index, grouping_source: "manual",
          name_source: "user", merged_into_id: null, flag: null, exposures: moved,
        };
        const srcIdx = ld.samples.indexOf(src);
        ld.samples.splice(srcIdx + 1, 0, created);
        break;
      }
      return loads;
    }),
  request: (p) => api.splitSample(p.sampleId, p.exposureIds, p.name, buildAuthOpts(p)),
  // Invalidate-only: the real new sample id arrives only via the refetch; the
  // own-op path never calls applyRemoteToCache, so we must invalidate here or
  // the negative-id placeholder is never reconciled (404 on a follow-up op).
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// dismiss_grouping_flag  (single-entity, DURABLE → grouping_flag_dismissed)
// sampleId in Input, NOT scope — one hook per experiment, not per row.
// ---------------------------------------------------------------------------
type DismissFlagInput = { sampleId: number; flagKind: "merge" | "split"; mergeWithSampleId?: number };
type DismissFlagScope = BaseScope;

export const dismissGroupingFlagMutator: Mutator<DismissFlagInput, DismissFlagScope, void> = {
  kind: "dismiss_grouping_flag",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s) { s.flag = null; break; } // optimistic: clear the suggestion
      }
      return loads;
    }),
  request: (p) => api.dismissGroupingFlag(
    p.sampleId,
    p.mergeWithSampleId !== undefined
      ? { flag_kind: p.flagKind, merge_with_sample_id: p.mergeWithSampleId }
      : { flag_kind: p.flagKind },
    buildAuthOpts(p),
  ),
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// undo_dismiss_grouping_flag  (DURABLE inverse of dismiss; re-shows the flag)
// Symmetric with dismissGroupingFlagMutator: suppress ↔ re-show.
// The optimistic patch restores the flag from the input (the page supplies the
// original flag snapshot from the server roll-up before dismiss fired); on
// success we invalidate loads to pick up the server-authoritative flag state.
// ---------------------------------------------------------------------------
export interface UndoDismissFlagInput {
  sampleId: number;
  /** The original flag to restore optimistically while the request is in-flight. */
  originalFlag: LoadSample["flag"];
}
type UndoDismissFlagScope = BaseScope;

export const undoDismissGroupingFlagMutator: Mutator<UndoDismissFlagInput, UndoDismissFlagScope, void> = {
  kind: "undo_dismiss_grouping_flag",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s) { s.flag = p.originalFlag; break; } // optimistic: re-show the flag
      }
      return loads;
    }),
  request: (p) => api.undoDismissGroupingFlag(p.sampleId, buildAuthOpts(p)),
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};
