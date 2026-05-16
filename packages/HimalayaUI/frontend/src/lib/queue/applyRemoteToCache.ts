import type { QueryClient } from "@tanstack/react-query";
import type { SseEvent, CurationPostState } from "./types";
import type {
  Peak, GroupEntry, Exposure, Sample, SampleMessage,
  ComparisonMessage, ComparisonSummary, Comparison,
} from "../../api";
import { queryKeys } from "../../queries";
import { peakQTol } from "./peakQTol";

/**
 * Write `post_state.indices` and `post_state.analysis_inputs_hash` into the
 * local cache without applying any per-kind body. Used by the foreign-event
 * path (via `applyRemoteToCache`) AND by the own-op confirmation paths in
 * `replayCoordinator` — own ops need post_state propagation too, otherwise
 * the indices cache stays frozen at the pre-mutation `inputs_hash` and the
 * StaleIndicesBanner sticks until a hard refetch.
 */
export function applyPostStateOnly(remote: SseEvent, qc: QueryClient): void {
  if (!remote.post_state) return;
  // `post_state` is a kind-keyed union (see SseEvent in types.ts): curation
  // frames (peak_*, analyze_run) carry a CurationPostState; comparison frames
  // (comparison_created/_submitted) carry a full Comparison. This function is
  // called unconditionally for ALL kinds from replayCoordinator's own-tab
  // paths (Case 1 deferred-match + self-echo guard), so an own-tab comparison
  // submit DOES reach here. It only writes the curation caches, so bail
  // unless the post_state actually is a CurationPostState — detected by its
  // `indices` array. Without this guard the `as CurationPostState` cast reads
  // `indices`/`analysis_inputs_hash` as `undefined` off a Comparison and
  // clobbers `queryKeys.exposure(entity_id)`; comparison ids and exposure ids
  // share the integer namespace, so a colliding cached exposure row would
  // have its `analysis_inputs_hash` wiped (falsely tripping StaleIndicesBanner
  // or breaking a later mutation's expected-hash check).
  const ps = remote.post_state as CurationPostState;
  if (!Array.isArray(ps.indices)) return;
  const id = remote.entity_id;
  qc.setQueryData(queryKeys.indices(id), ps.indices);
  qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
    old
      ? { ...old, analysis_inputs_hash: ps.analysis_inputs_hash }
      : old);
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
    case "index_confirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      const cached = qc.getQueryData<GroupEntry[]>(queryKeys.groups(id));
      // Issue #37 Bug 1c: the foreign tab may have just triggered
      // `ensure_custom_group!`, which mints a new custom group id and
      // demotes the auto group. Local tabs won't have that group cached;
      // the surgical splice would silently miss the foreign confirmation.
      // Invalidate to refetch the canonical group structure.
      if (!cached || !cached.some((g) => g.id === groupId)) {
        qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
        break;
      }
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), cached.map((g) =>
        g.id === groupId ? { ...g, members: [...g.members, indexId] } : g));
      break;
    }
    case "index_unconfirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      const cached = qc.getQueryData<GroupEntry[]>(queryKeys.groups(id));
      if (!cached || !cached.some((g) => g.id === groupId)) {
        qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
        break;
      }
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), cached.map((g) =>
        g.id === groupId
          ? { ...g, members: g.members.filter((m) => m !== indexId) }
          : g));
      break;
    }
    case "speculative_created":
    case "speculative_deleted": {
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
      break;
    }
    case "set_exposure_status": {
      qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
        old ? { ...old, status: payload?.status as Exposure["status"] } : old);
      break;
    }
    case "post_message": {
      // The same OpKind covers two distinct chat threads. The route handler
      // tags `entity_type` so the merger can pick the right cache key:
      //   - `sample_message`     → `routes_messages.jl` (sample chat)
      //   - `comparison_message` → `routes_comparisons.jl` (comparison chat)
      // entity_id on the SSE frame is the message id; the parent id rides in
      // the payload.
      if (remote.entity_type === "comparison_message") {
        const comparisonId = payload?.comparison_id as number;
        qc.setQueryData<ComparisonMessage[]>(
          queryKeys.comparisonMessages(comparisonId),
          (old = []) => {
            const incoming = payload as unknown as ComparisonMessage;
            if (old.some((m) => m.id === incoming.id)) return old;
            return [...old, incoming];
          });
      } else {
        // Default + explicit `sample_message`: legacy behavior.
        const sampleId = payload?.sample_id as number;
        qc.setQueryData<SampleMessage[]>(queryKeys.messages(sampleId), (old = []) => {
          // Dedupe in case the same SSE arrives twice (own-op late echo races
          // post-onSuccess writes); message ids are server-assigned positives.
          const incoming = payload as unknown as SampleMessage;
          if (old.some((m) => m.id === incoming.id)) return old;
          return [...old, incoming];
        });
      }
      break;
    }
    case "comparison_created":
    case "comparison_submitted": {
      // Foreign tab created or re-submitted a comparison. As of Compare UX
      // A-5 Step 5b the SSE frame carries `post_state =
      // fetch_comparison_with_members(db, id)` — the exact same projection
      // `GET /api/comparisons/:id` returns, including server-computed
      // `forked_from_title` and the persisted view_* fields. Splicing it
      // straight into the cache is now both safe (no client-side
      // re-derivation to drift) AND required (A-9: view_* must land without
      // a refetch round-trip so the detail page reflects the author's view).
      if (remote.post_state != null) {
        const post = remote.post_state as Comparison;
        qc.setQueryData(queryKeys.comparison(id), post);
        if (Array.isArray(post.members)) {
          qc.setQueryData(queryKeys.comparisonMembers(id), post.members);
        }
      } else {
        // Fallback: a pre-A-5 frame without post_state. Keep the
        // invalidate safety net so a partial deploy can't strand the cache.
        qc.invalidateQueries({ queryKey: queryKeys.comparison(id) });
        qc.invalidateQueries({ queryKey: queryKeys.comparisonMembers(id) });
      }
      // The listing cache always invalidates — too many denormalised
      // projection fields (member_count, member_phases, last_event_at,
      // has_stale_members) for a manual splice, and membership-derived
      // listings change in either direction on any submit.
      qc.invalidateQueries({ queryKey: ["comparisons"] });
      break;
    }
    case "comparison_deleted": {
      // Remove the entity caches — refetching a deleted resource 404s and
      // leaves stale `isError: true` state. Spec line 345.
      qc.removeQueries({ queryKey: queryKeys.comparison(id) });
      qc.removeQueries({ queryKey: queryKeys.comparisonMembers(id) });
      qc.removeQueries({ queryKey: queryKeys.comparisonMessages(id) });
      qc.removeQueries({ queryKey: queryKeys.comparisonForks(id) });
      // Filter the id out of every cached listing (prefix walk).
      qc.setQueriesData<ComparisonSummary[]>(
        { queryKey: ["comparisons"] },
        (old) => (old ? old.filter((c) => c.id !== id) : old),
      );
      break;
    }
    case "comparison_pinned":
    case "comparison_unpinned": {
      // Pin/unpin SSE drives cross-tab fanout. The frontend's
      // `comparisonPins` cache is global per-tab (only ever stores the
      // current user's pin set), and the SSE self-echo filter via
      // `client_id` already discards events the originating tab emitted.
      // What's left are foreign-tab pins from the same user → invalidate
      // so the next read gets the canonical list. Pins from OTHER users
      // also fan out (every connected client sees every event), but their
      // user_id mismatches and the route's per-user filter excludes them
      // on read. Invalidating in that case is a harmless no-op refetch.
      qc.invalidateQueries({ queryKey: queryKeys.comparisonPins });
      break;
    }
    case "add_tag":
    case "remove_tag": {
      const parentKey = remote.entity_type === "sample"
        ? queryKeys.samples(payload?.experiment_id as number)
        : queryKeys.exposures(payload?.sample_id as number);
      qc.invalidateQueries({ queryKey: parentKey });
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
      break;
    }
    case "analyze_run": {
      applyPostState();
      break;
    }
    default: {
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
    }
  }
}
