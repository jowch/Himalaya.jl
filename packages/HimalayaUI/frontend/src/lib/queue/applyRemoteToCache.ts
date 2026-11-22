import { refreshPeaks } from "./refreshPeaks";
import type { QueryClient } from "@tanstack/react-query";
import type { SseEvent, CurationPostState, AssignmentPostState, SampleCallPostState } from "./types";
import { writeSampleCall, invalidateOtherFrameAssignments, invalidateSampleCalls } from "./sampleCall";
import type {
  Peak, Exposure, Sample,
  Assignment, Experiment,
} from "../../api";
import { queryKeys } from "../../queries";
import { mergeExperimentPatch } from "./mutators/experiment";
import { mergeSetupPatch, type SetupFields } from "./mutators/setup";
import { peakQTol } from "./peakQTol";
import { announceAssignmentDropped } from "./assignmentDropped";
import { INGEST_STAGES } from "../ingestStages";

/** The one stage that predates the persist txn — see invalidateIngestFrameCache. */
const INGEST_DISCOVERY_STAGE = INGEST_STAGES[0];

/**
 * Write the `{state, members}` assignment envelope into the assignment cache.
 * Single writer shared by the assignment_* branch in `applyRemoteToCache`
 * and the curation-frame path in `applyPostStateOnly` (F-WIPE W2: peak_*
 * frames now carry the same envelope, serialized by the same backend
 * helper) — one place to keep the cache row shape honest.
 */
function writeAssignmentFromPostState(
  exposureId: number,
  assignment: Pick<Assignment, "state" | "members">,
  qc: QueryClient,
): void {
  qc.setQueryData<Assignment>(queryKeys.assignment(exposureId), {
    exposure_id: exposureId,
    state: assignment.state,
    members: assignment.members,
  });
}

export function applyPostStateOnly(remote: SseEvent, qc: QueryClient): void {
  if (!remote.post_state) return;
  const ps = remote.post_state as CurationPostState;
  if (!Array.isArray(ps.indices)) return;
  const id = remote.entity_id;
  qc.setQueryData(queryKeys.indices(id), ps.indices);
  // F-WIPE W2: every reanalyzing frame (peak_* AND analyze_run) carries the
  // re-attached assignment inside the curation post_state (same {state,
  // members} envelope as the assignment_* frames), plus assignment_dropped
  // when non-empty. Splice it; when ABSENT (pre-W1 backend) do NOT
  // invalidate or touch the assignment cache, preserving the previous
  // behavior.
  if (ps.assignment) {
    writeAssignmentFromPostState(id, ps.assignment, qc);
  }
  // Announce dropped call members. Consequential state change → visible
  // toast; the Toast surface is itself the aria-live region (role="alert"
  // for warnings), so no separate announce() — see assignmentDropped.ts.
  // Exactly-once: see the path enumeration in this function's docstring.
  announceAssignmentDropped(ps.assignment_dropped);
}

/**
 * Shared cache-invalidation helper for ingest-frame side-effects. Exported so
 * `App.tsx`'s SSE listener can call the same invalidations without duplicating
 * them inline (the listener still owns the Zustand store write; this helper is
 * pure cache-only). `isComplete=true` also refetches the experiment detail row
 * (so `ingest_status` transitions from "analyzing" → "complete").
 *
 * `stage` skips the loads/samples refetch during `"discovery"`. Nothing is
 * committed then — `scan_and_group!` inserts every load/sample/exposure in ONE
 * transaction that has not opened yet — so the refetch cannot observe anything
 * new, once per progress frame, against the slowest payload on the page. (On a
 * RESCAN the rows do already exist, so this is "returns what it already had",
 * not "returns empty" — either way the refetch is wasted until the txn commits.) From `"analyzing"` onward the txn HAS committed and
 * `analyze_exposure!` mutates per-exposure data, so the refetch is real work.
 *
 * An ABSENT stage keeps the original always-invalidate behavior — an
 * `ingest_started` frame, or a backend predating stage reporting, must not
 * silently lose its refresh.
 */
export function invalidateIngestFrameCache(
  qc: QueryClient,
  expId: number,
  isComplete: boolean,
  stage?: string,
): void {
  // isComplete wins over stage: a terminal frame must always refresh, whatever
  // stage it happens to carry. Without this the two callers disagreed --
  // App.tsx forwards `stage` for all four ingest kinds, this file's own
  // ingest_complete arm passes none.
  if (isComplete || stage !== INGEST_DISCOVERY_STAGE) {
    qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
    qc.invalidateQueries({ queryKey: queryKeys.experimentCorpusSamples(expId) });
  }
  if (isComplete) {
    qc.invalidateQueries({ queryKey: queryKeys.experiment(expId) });
    qc.invalidateQueries({ queryKey: queryKeys.problems(expId) });
    qc.invalidateQueries({ queryKey: queryKeys.labels(expId) });
  }
}

/**
 * Apply a remote SSE event to the local query cache. Per-kind logic mirrors
 * the spec's "replay-without-refetch where post_state covers it; refetch
 * fallback where the event payload is insufficient or update is rare."
 *
 * Split out of replayCoordinator.ts because this dispatcher will grow as M2.1+
 * activates more event kinds; replayCoordinator.ts stays focused on queue
 * orchestration (deferred resolution + rollback/rerun).
 */
export function applyRemoteToCache(remote: SseEvent, qc: QueryClient): void {
  const id = remote.entity_id;
  const payload = remote.payload as Record<string, unknown> | undefined;

  const applyPostState = (): void => applyPostStateOnly(remote, qc);

  switch (remote.kind) {
    case "peak_added": {
      // Server includes peak_curation_id (= the real DB id) in the payload —
      // use it as the row's id so a foreign-tab insert can later be deleted /
      // patched without a 404. (Pre-fix: used -remote.id, the EVENT id, which
      // didn't correspond to any real peak row — issue #2 from PR review.)
      const peakId = payload?.peak_curation_id as number | undefined;
      if (peakId === undefined) {
        // Defensive fallback: if a future server change drops peak_curation_id,
        // invalidate so the next read replaces the placeholder with the real
        // row rather than leaving a phantom in the cache.
        qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
        applyPostState();
        break;
      }
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) => {
        // Idempotent insert: dedupe against an existing row (own-tab SSE echo
        // arriving after onSuccess already wrote the canonical row).
        // Source-scoped: auto and curation tables have independent id sequences
        // and can collide on the wire. Same invariant as peakAdd mutator.
        if (old.some((p) => p.id === peakId && p.source === "manual")) return old;
        return [
          ...old,
          {
            id: peakId,
            exposure_id: id,
            q: payload?.q as number,
            intensity: null,
            prominence: null,
            sharpness: null,
            source: "manual",
            excluded: false,
          },
        ];
      });
      refreshPeaks(qc, id);
      applyPostState();
      break;
    }
    case "peak_excluded": {
      // Prefer auto_peak_id when present — id match is unambiguous, q match
      // can mis-pick when two auto peaks are within tolerance of each other
      // (suggestion #1 from PR review).
      const autoPeakId = payload?.auto_peak_id as number | undefined;
      const targetQ = payload?.q as number;
      const tol = peakQTol(targetQ);
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) => {
          // Only auto peaks can be excluded; scope by source to avoid flipping
          // a manual peak that happens to share an id with the targeted auto peak.
          const matches = autoPeakId !== undefined
            ? p.id === autoPeakId && p.source === "auto"
            : p.source === "auto" && Math.abs(p.q - targetQ) < tol;
          return matches ? { ...p, excluded: true } : p;
        }));
      refreshPeaks(qc, id);
      applyPostState();
      break;
    }
    case "peak_unexcluded": {
      const autoPeakId = payload?.auto_peak_id as number | undefined;
      const targetQ = payload?.q as number;
      const tol = peakQTol(targetQ);
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) => {
          const matches = autoPeakId !== undefined
            ? p.id === autoPeakId && p.source === "auto"
            : p.source === "auto" && Math.abs(p.q - targetQ) < tol;
          return matches ? { ...p, excluded: false } : p;
        }));
      refreshPeaks(qc, id);
      applyPostState();
      break;
    }
    case "peak_removed": {
      // Filter by peak_curation_id when present; q-tolerance is the fallback
      // for events emitted before the id was added to the payload (issue #1).
      const removedId = payload?.peak_curation_id as number | undefined;
      const targetQ = payload?.q as number | undefined;
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) => {
        if (removedId !== undefined) {
          // Only manual peaks are removable; scope by source so a colliding
          // auto peak id can't be dropped by a foreign-tab remove event.
          return old.filter((p) => !(p.id === removedId && p.source === "manual"));
        }
        if (targetQ !== undefined) {
          const tol = peakQTol(targetQ);
          return old.filter((p) => !(p.source === "manual" && Math.abs(p.q - targetQ) < tol));
        }
        // Neither id nor q — payload is unusable; refetch.
        qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
        return old;
      });
      refreshPeaks(qc, id);
      applyPostState();
      break;
    }
    case "sample_call_set_state": {
      const ps = remote.post_state as SampleCallPostState | undefined;
      if (ps?.call) writeSampleCall(qc, id, ps.call);
      else qc.invalidateQueries({ predicate: (q) => q.queryKey[0] === "exposure" && q.queryKey[2] === "assignment" });
      invalidateSampleCalls(qc);
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
    case "assignment_add":
    case "assignment_remove":
    case "assignment_set_state": {
      // Assignment frames carry a DISTINCT post_state — {assignment:{state,
      // members}} with no top-level `indices` key (so applyPostStateOnly bails
      // for these kinds, never touching the curation caches). Patch the
      // assignment cache directly from post_state; invalidate as the fallback
      // for a (pre-Plan-D) frame that lacks the envelope.
      const ps = remote.post_state as AssignmentPostState | undefined;
      if (ps?.assignment) {
        writeAssignmentFromPostState(id, ps.assignment, qc);
      } else {
        qc.invalidateQueries({ queryKey: queryKeys.assignment(id) });
      }
      if (remote.kind === "assignment_add") invalidateOtherFrameAssignments(qc, id);
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      invalidateSampleCalls(qc);
      break;
    }
    case "speculative_created": {
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      break;
    }
    case "speculative_deleted": {
      // Also the assignment: `assignment_members.index_id` is ON DELETE
      // CASCADE (db.jl:223), so the route's DELETE silently removes the index
      // from the durable call. `speculative_deleted` carries no post_state, so
      // a foreign tab that only invalidated `indices` would keep a member id
      // pointing at an index that no longer exists — a phantom cart block.
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.assignment(id) });
      break;
    }
    case "set_exposure_status": {
      const status = payload?.status as Exposure["status"];
      qc.setQueriesData<Exposure[]>(
        { predicate: (q) => q.queryKey[0] === "sample" && q.queryKey[2] === "exposures" },
        (list) => list?.some((e) => e.id === id)
          ? list.map((e) => (e.id === id ? { ...e, status } : e))
          : list,
      );
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
    case "update_sample": {
      qc.setQueriesData<Sample[]>(
        { predicate: (q) => q.queryKey[0] === "experiment" && q.queryKey[2] === "samples" },
        (list) => list?.some((s) => s.id === id)
          ? list.map((s) => (s.id === id ? { ...s, ...(payload ?? {}) } : s))
          : list,
      );
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
    case "analyze_run": {
      applyPostState();
      if (payload?.trace_hash_before !== payload?.trace_hash_after) {
        qc.invalidateQueries({ queryKey: queryKeys.trace(id) });
        qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      }
      break;
    }
    case "ingest_started":
    case "ingest_progress": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) {
        invalidateIngestFrameCache(qc, expId, false, payload?.stage as string | undefined);
      }
      break;
    }
    case "ingest_complete":
    case "ingest_failed": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) invalidateIngestFrameCache(qc, expId, true);
      break;
    }
    case "sample_bridged":
    case "sample_unbridged": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.bridges(expId) });
      invalidateSampleCalls(qc);
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
    case "compares_changed": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.compares(expId) });
      break;
    }
    case "tag_color_set":
      qc.invalidateQueries({ queryKey: queryKeys.labelColors() });
      break;
    case "experiment_updated":
      mergeExperimentPatch(qc, id, (payload ?? {}) as Partial<Experiment>);
      qc.invalidateQueries({ queryKey: queryKeys.experiment(id) });
      qc.invalidateQueries({ queryKey: queryKeys.experiments });
      break;
    case "setup_updated": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId === undefined) break;
      const fields: SetupFields = {};
      if (typeof payload?.name === "string") fields.name = payload.name;
      if (payload !== undefined && "display_name" in payload) fields.display_name = (payload.display_name as string | null) ?? null;
      if (typeof payload?.included === "boolean") fields.included = payload.included;
      mergeSetupPatch(qc, expId, id, fields);
      qc.invalidateQueries({ queryKey: queryKeys.experiment(expId) });
      qc.invalidateQueries({ queryKey: queryKeys.experiments });
      break;
    }
    case "label_added":
    case "label_removed":
    case "label_suggestion_restored": {
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.labels(expId) });
      break;
    }
    default: {
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
    }
  }
}
