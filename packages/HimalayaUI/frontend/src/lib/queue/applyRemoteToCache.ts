import type { QueryClient } from "@tanstack/react-query";
import type { SseEvent, CurationPostState, AssignmentPostState } from "./types";
import type {
  Peak, Exposure, Sample, SampleMessage,
  Series, Assignment,
} from "../../api";
import { queryKeys } from "../../queries";
import { peakQTol } from "./peakQTol";
import { announceAssignmentDropped } from "./assignmentDropped";

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

/**
 * Write `post_state.indices` and `post_state.analysis_inputs_hash` into the
 * local cache without applying any per-kind body. Used by the foreign-event
 * path (via `applyRemoteToCache`) AND by the own-op confirmation paths in
 * `replayCoordinator` — own ops need post_state propagation too, otherwise
 * the indices cache stays frozen at the pre-mutation `inputs_hash` and the
 * stale-indices alert sticks until a hard refetch.
 *
 * F-WIPE W2: this is ALSO where the W1 assignment envelope on peak_* frames
 * is consumed, deliberately. Every curation frame is processed by exactly
 * one of three mutually-exclusive paths per tab, and ALL of them funnel
 * through this function:
 *
 *   1. Own-op SSE-wins (replayCoordinator Case 1, deferred match) — calls
 *      applyPostStateOnly directly, then resolves the deferred (the HTTP
 *      request is aborted; its response is never processed).
 *   2. Own-op HTTP-wins (self-echo guard: client_id matches, deferred
 *      already cleared) — calls applyPostStateOnly directly. The mutator's
 *      onSuccess (HTTP path) never reads post_state, so nothing double-fires.
 *   3. Foreign frame — applyRemoteToCache's peak_* and analyze_run branches
 *      call applyPostState() → here.
 *
 * Placing the assignment write AND the assignment_dropped announcement here
 * therefore guarantees exactly-once semantics for both the editing tab and
 * foreign tabs, with no per-branch duplication.
 */
export function applyPostStateOnly(remote: SseEvent, qc: QueryClient): void {
  if (!remote.post_state) return;
  // `post_state` is a kind-keyed union (see SseEvent in types.ts): curation
  // frames (peak_*, analyze_run) carry a CurationPostState; series-commit
  // frames (series_plate_committed) carry a full Series. This function is
  // called unconditionally for ALL kinds from replayCoordinator's own-tab
  // paths (Case 1 deferred-match + self-echo guard), so an own-tab series
  // commit DOES reach here. It only writes the curation caches, so bail
  // unless the post_state actually is a CurationPostState — detected by its
  // `indices` array. Without this guard the `as CurationPostState` cast reads
  // `indices`/`analysis_inputs_hash` as `undefined` off a Series and
  // clobbers `queryKeys.exposure(entity_id)`; series ids and exposure ids
  // share the integer namespace, so a colliding cached exposure row would
  // have its `analysis_inputs_hash` wiped (falsely tripping the stale-indices alert
  // or breaking a later mutation's expected-hash check).
  const ps = remote.post_state as CurationPostState;
  if (!Array.isArray(ps.indices)) return;
  const id = remote.entity_id;
  qc.setQueryData(queryKeys.indices(id), ps.indices);
  qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
    old
      ? { ...old, analysis_inputs_hash: ps.analysis_inputs_hash }
      : old);
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
 */
export function invalidateIngestFrameCache(
  qc: QueryClient,
  expId: number,
  isComplete: boolean,
): void {
  qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
  qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
  if (isComplete) qc.invalidateQueries({ queryKey: queryKeys.experiment(expId) });
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
      applyPostState();
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
      qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
        old ? { ...old, status: payload?.status as Exposure["status"] } : old);
      break;
    }
    case "post_message": {
      // Sample chat only (the comparison-message thread was retired with the
      // Compare page). entity_id on the SSE frame is the message id; the
      // parent sample id rides in the payload.
      const sampleId = payload?.sample_id as number;
      qc.setQueryData<SampleMessage[]>(queryKeys.messages(sampleId), (old = []) => {
        // Dedupe in case the same SSE arrives twice (own-op late echo races
        // post-onSuccess writes); message ids are server-assigned positives.
        const incoming = payload as unknown as SampleMessage;
        if (old.some((m) => m.id === incoming.id)) return old;
        return [...old, incoming];
      });
      break;
    }
    case "series_created":
    case "series_recipe_updated": {
      // Neither kind carries post_state (master plan §5.2); the SSE payload's
      // series_samples entries are id-less and series_samples.id is
      // replay-volatile (§11), so there is no safe surgical splice.
      // Invalidate-only — the next read refetches the canonical projection.
      qc.invalidateQueries({ queryKey: queryKeys.series(id) });
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
      break;
    }
    case "series_deleted": {
      // Remove the detail cache — refetching a deleted resource 404s and
      // leaves stale `isError` state. Filter the id out of the listing.
      qc.removeQueries({ queryKey: queryKeys.series(id) });
      qc.setQueriesData<{ id: number }[]>(
        { queryKey: queryKeys.seriesList },
        (old) => (old ? old.filter((s) => s.id !== id) : old),
      );
      break;
    }
    case "series_plate_committed": {
      // The one series event carrying a post_state envelope (master plan
      // §5.2): post_state IS the fetch_series_with_plate projection. Splice it
      // straight into the detail cache; invalidate the listing (denormalised
      // member_count / has_stale_members / last_event_at fields).
      if (remote.post_state != null) {
        qc.setQueryData(queryKeys.series(id), remote.post_state as Series);
      } else {
        qc.invalidateQueries({ queryKey: queryKeys.series(id) });
      }
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
      break;
    }
    case "add_tag":
    case "remove_tag":
    case "edit_tag": {
      if (remote.entity_type === "sample") {
        // A sample tag appears in two cached projections: the per-experiment
        // samples list AND the corpus-wide contact-sheet list. The route
        // always emits experiment_id, so the experiment key still invalidates
        // (a harmless no-op if not cached); the corpusSamples invalidation is
        // what refreshes the contact sheet from a foreign tab (#159).
        qc.invalidateQueries({
          queryKey: queryKeys.samples(payload?.experiment_id as number),
        });
        qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
        // I3.4 (#174): the scoping surface (/series/new) reads two more corpus
        // projections off the same tags — the distinct (key,value) proposal
        // source and the picker projection (which carries per-sample tags). A
        // foreign scoping write must refresh both so a peer's proposal isn't
        // stale (master plan §11 — add_tag fan-out re-audit). Tri-scope-safe:
        // this is inside the `entity_type === "sample"` arm, so exposure/corpus
        // tag routing is untouched.
        qc.invalidateQueries({ queryKey: queryKeys.corpusSampleTags });
        qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
      } else {
        qc.invalidateQueries({
          queryKey: queryKeys.exposures(payload?.sample_id as number),
        });
      }
      break;
    }
    case "update_sample": {
      qc.setQueryData(queryKeys.sample(id), (old: Sample | undefined) =>
        old ? { ...old, ...(payload ?? {}) } : old);
      break;
    }
    case "select_exposure": {
      const sampleId = payload?.sample_id as number;
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) });
      // The picker's indexing_exposure_id derives from exposure selection
      // (highest-id selected, else highest-id), and the builder's Confirm
      // resolves its plate through that projection (BU-RECIPENOOP) - a stale
      // picker after a foreign re-selection would silently commit the
      // PREVIOUS representative exposure.
      qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
      break;
    }
    case "analyze_run": {
      applyPostState();
      break;
    }
    case "ingest_started":
    case "ingest_progress":
    case "ingest_failed": {
      // Broadcast-only progress (spec §9.3): rides the curation channel, carries
      // a positive experiment_id in the payload (NEVER a sentinel entity_id).
      // Read experiment_id from the payload, not remote.entity_id. Cache effect
      // is invalidation-only (the ingestInFlight store write lives in the
      // separate App.tsx listener — applyRemoteToCache stays pure).
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) invalidateIngestFrameCache(qc, expId, false);
      break;
    }
    case "ingest_complete": {
      // Authoritative terminal frame (the 64-slot channel may drop progress
      // frames at 680-exposure scale; treat complete as the source of truth).
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) invalidateIngestFrameCache(qc, expId, true);
      break;
    }
    case "sample_renamed": {
      // Single-entity, payload-derivable: splice the renamed label into the
      // sample cache + refresh the corpus/experiment listings. entity_id is the
      // sample id; the new label rides payload.name.
      const newName = payload?.name as string | undefined;
      if (newName !== undefined) {
        qc.setQueryData<Sample>(queryKeys.sample(id), (old) =>
          old ? { ...old, name: newName } : old);
      }
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
      break;
    }
    case "exposure_moved": {
      // One exposure left from_sample_id for sample_id (the destination).
      // Invalidate both sides' exposure lists + the loads roll-up (which
      // re-derives the grouping view).
      const dest = payload?.sample_id as number | undefined;
      const from = payload?.from_sample_id as number | undefined;
      if (from !== undefined) qc.invalidateQueries({ queryKey: queryKeys.exposures(from) });
      if (dest !== undefined) qc.invalidateQueries({ queryKey: queryKeys.exposures(dest) });
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
      break;
    }
    case "sample_created":
    case "sample_split":
    case "grouping_flag_dismissed": {
      // Orchestrations that create sample rows / change grouping whose new ids
      // aren't worth a surgical splice (the series_created precedent). Refetch
      // the loads roll-up + both sample listings. MUST be before `default:` so
      // entity_id (a sample id) never poisons peaks(id)/indices(id). There is
      // NO sample_merged event — merging is grouping_flag_dismissed + the move.
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) {
        qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
        qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
      }
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
    default: {
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
    }
  }
}
