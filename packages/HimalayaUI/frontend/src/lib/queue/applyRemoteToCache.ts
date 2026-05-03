import type { QueryClient } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import type { Peak, GroupEntry, Exposure, Sample, SampleMessage } from "../../api";
import { queryKeys } from "../../queries";

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
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) => [
        ...old,
        // Negative id placeholder per the optimistic-id invariant in types.ts.
        // The remote event represents a manually-added peak; intensity/
        // prominence/sharpness are null because the server hasn't run analysis
        // yet. Peak's interface allows null for these fields.
        {
          id: -remote.id,
          exposure_id: id,
          q: payload?.q as number,
          intensity: null,
          prominence: null,
          sharpness: null,
          source: "manual",
          excluded: false,
        },
      ]);
      applyPostState();
      break;
    }
    case "peak_excluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) =>
          Math.abs(p.q - (payload?.q as number)) < 1e-6
            ? { ...p, excluded: true }
            : p));
      applyPostState();
      break;
    }
    case "peak_unexcluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) =>
          Math.abs(p.q - (payload?.q as number)) < 1e-6
            ? { ...p, excluded: false }
            : p));
      applyPostState();
      break;
    }
    case "peak_removed": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.filter((p) => Math.abs(p.q - (payload?.q as number)) >= 1e-6));
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
      // Forward-scaffolded: route still uses log_action!.
      qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
        old ? { ...old, status: payload?.status as Exposure["status"] } : old);
      break;
    }
    case "post_message": {
      // Forward-scaffolded: route still uses log_action!.
      const sampleId = payload?.sample_id as number;
      qc.setQueryData<SampleMessage[]>(queryKeys.messages(sampleId), (old = []) => [
        ...old,
        payload as unknown as SampleMessage,
      ]);
      break;
    }
    case "add_tag":
    case "remove_tag": {
      // Forward-scaffolded: route still uses log_action!.
      const parentKey = remote.entity_type === "sample"
        ? queryKeys.samples(payload?.experiment_id as number)
        : queryKeys.exposures(payload?.sample_id as number);
      qc.invalidateQueries({ queryKey: parentKey });
      break;
    }
    case "update_sample": {
      // Forward-scaffolded: route still uses log_action!.
      qc.setQueryData(queryKeys.sample(id), (old: Sample | undefined) =>
        old ? { ...old, ...(payload ?? {}) } : old);
      break;
    }
    case "select_exposure": {
      // Forward-scaffolded: route still uses log_action!.
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
