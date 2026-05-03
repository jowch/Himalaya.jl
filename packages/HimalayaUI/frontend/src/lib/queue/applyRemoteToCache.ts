import type { QueryClient } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import type { Peak, GroupEntry, Exposure, Sample, SampleMessage } from "../../api";
import { queryKeys } from "../../queries";
import { peakQTol } from "./peakQTol";

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

  const applyPostState = (): void => {
    if (!remote.post_state) return;
    qc.setQueryData(queryKeys.indices(id), remote.post_state.indices);
    qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
      old
        ? { ...old, analysis_inputs_hash: remote.post_state!.analysis_inputs_hash }
        : old);
  };

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
        if (old.some((p) => p.id === peakId)) return old;
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
          const matches = autoPeakId !== undefined
            ? p.id === autoPeakId
            : Math.abs(p.q - targetQ) < tol;
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
            ? p.id === autoPeakId
            : Math.abs(p.q - targetQ) < tol;
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
          return old.filter((p) => p.id !== removedId);
        }
        if (targetQ !== undefined) {
          const tol = peakQTol(targetQ);
          return old.filter((p) => Math.abs(p.q - targetQ) >= tol);
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
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), (old = []) =>
        old.map((g) =>
          g.id === groupId ? { ...g, members: [...g.members, indexId] } : g));
      break;
    }
    case "index_unconfirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), (old = []) =>
        old.map((g) =>
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
      // entity_id on the SSE frame is the sample_message id; the message's
      // sample_id rides in the payload (set by the route handler).
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
