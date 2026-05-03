/**
 * Trivial-mutation mutators (M2.1) — simple field-level optimistic updates
 * with no peak-set / index-set side effects. Each mutator wires:
 * - `onMutate`: writes the optimistic effect, returns a `restore` closure that
 *   reverts the cache to its pre-mutation snapshot on rollback or HTTP error
 * - `request`: thin wrapper around the corresponding `api.*` call
 * - `onSuccess`: replaces the optimistic placeholder (where applicable) with
 *   the server's authoritative response. SSE-driven cross-tab sync is handled
 *   separately by `applyRemoteToCache`; we deliberately avoid `invalidateQueries`
 *   here so a mutating tab does not pay a refetch round-trip on every action.
 */
import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type {
  Sample, SampleTag, Exposure, ExposureTag, SampleMessage, AuthOpts,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, OpPayload, RollbackContext } from "../types";

// ---------------------------------------------------------------------------
// `useQueueMutation` flat-spreads the scope object and the per-call input into
// the payload it constructs (`{ kind, clientOpId, ...scope, ...input }`),
// even though its type signature describes the payload as `OpPayload<TInput>`.
// We declare each mutator's payload as `OpPayload<TInput>` to satisfy the
// framework's constraint, then narrow to a per-mutator `Flat<…>` type via a
// single cast at the top of each callback to access the merged scope fields.
// ---------------------------------------------------------------------------

interface BaseScope {
  username: string | undefined;
  clientId: string;
}

type Flat<TInput, TScope extends object> = OpPayload<TInput> & TScope & TInput & BaseScope;

const flat = <TInput, TScope extends object>(
  p: OpPayload<TInput>,
): Flat<TInput, TScope> => p as unknown as Flat<TInput, TScope>;

function buildAuthOpts(p: {
  username?: string | undefined;
  clientId?: string | undefined;
  clientOpId?: string | undefined;
}): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/** Generate a unique negative id for an optimistic placeholder row. */
let optimisticCounter = 0;
function nextOptimisticId(): number {
  optimisticCounter += 1;
  return -(Date.now() * 1000 + optimisticCounter);
}

// ---------------------------------------------------------------------------
// update_sample
// ---------------------------------------------------------------------------

export type UpdateSampleInput = { name?: string; notes?: string };
type UpdateSampleScope = { experimentId: number; sampleId: number };

export const updateSampleMutator: Mutator<OpPayload<UpdateSampleInput>, Sample> = {
  kind: "update_sample",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<UpdateSampleInput, UpdateSampleScope>(raw);
    const samplesKey = queryKeys.samples(p.experimentId);
    const sampleKey = queryKeys.sample(p.sampleId);
    const prevList = qc.getQueryData<Sample[]>(samplesKey);
    const prevSingle = qc.getQueryData<Sample>(sampleKey);
    const patch = patchOf(p);

    if (prevList) {
      qc.setQueryData<Sample[]>(samplesKey, prevList.map((s) =>
        s.id === p.sampleId ? { ...s, ...patch } : s,
      ));
    }
    if (prevSingle && prevSingle.id === p.sampleId) {
      qc.setQueryData<Sample>(sampleKey, { ...prevSingle, ...patch });
    }
    return {
      restore: () => {
        if (prevList !== undefined) qc.setQueryData(samplesKey, prevList);
        if (prevSingle !== undefined) qc.setQueryData(sampleKey, prevSingle);
      },
    };
  },
  request: (raw) => {
    const p = flat<UpdateSampleInput, UpdateSampleScope>(raw);
    return api.updateSample(p.sampleId, patchOf(p), buildAuthOpts(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flat<UpdateSampleInput, UpdateSampleScope>(raw);
    const samplesKey = queryKeys.samples(p.experimentId);
    const sampleKey = queryKeys.sample(p.sampleId);
    const list = qc.getQueryData<Sample[]>(samplesKey);
    if (list) {
      qc.setQueryData<Sample[]>(samplesKey, list.map((s) =>
        s.id === response.id ? response : s));
    }
    qc.setQueryData<Sample>(sampleKey, response);
  },
};

function patchOf(p: { name?: string; notes?: string }): { name?: string; notes?: string } {
  const out: { name?: string; notes?: string } = {};
  if (p.name !== undefined) out.name = p.name;
  if (p.notes !== undefined) out.notes = p.notes;
  return out;
}

// ---------------------------------------------------------------------------
// add_tag (sample)
// ---------------------------------------------------------------------------

export type AddSampleTagInput = { key: string; value: string };
type AddSampleTagScope = { experimentId: number; sampleId: number };

export const addSampleTagMutator: Mutator<OpPayload<AddSampleTagInput>, SampleTag> = {
  kind: "add_tag",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<AddSampleTagInput, AddSampleTagScope>(raw);
    const samplesKey = queryKeys.samples(p.experimentId);
    const prev = qc.getQueryData<Sample[]>(samplesKey);
    const placeholderId = nextOptimisticId();
    if (prev) {
      qc.setQueryData<Sample[]>(samplesKey, prev.map((s) =>
        s.id === p.sampleId
          ? { ...s, tags: [...s.tags, {
              id: placeholderId, key: p.key, value: p.value, source: "manual",
            }] }
          : s,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(samplesKey, prev);
      },
    };
  },
  request: (raw) => {
    const p = flat<AddSampleTagInput, AddSampleTagScope>(raw);
    return api.addSampleTag(p.sampleId, p.key, p.value, buildAuthOpts(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flat<AddSampleTagInput, AddSampleTagScope>(raw);
    const samplesKey = queryKeys.samples(p.experimentId);
    qc.setQueryData<Sample[]>(samplesKey, (list) => {
      if (!list) return list;
      return list.map((s) => {
        if (s.id !== p.sampleId) return s;
        const filtered = s.tags.filter((t) =>
          !(t.id < 0 && t.key === p.key && t.value === p.value)
          && t.id !== response.id,
        );
        return { ...s, tags: [...filtered, response] };
      });
    });
  },
};

// ---------------------------------------------------------------------------
// remove_tag (sample)
// ---------------------------------------------------------------------------

export type RemoveSampleTagInput = { tagId: number };
type RemoveSampleTagScope = { experimentId: number; sampleId: number };

export const removeSampleTagMutator: Mutator<OpPayload<RemoveSampleTagInput>, void> = {
  kind: "remove_tag",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<RemoveSampleTagInput, RemoveSampleTagScope>(raw);
    const samplesKey = queryKeys.samples(p.experimentId);
    const prev = qc.getQueryData<Sample[]>(samplesKey);
    if (prev) {
      qc.setQueryData<Sample[]>(samplesKey, prev.map((s) =>
        s.id === p.sampleId
          ? { ...s, tags: s.tags.filter((t) => t.id !== p.tagId) }
          : s,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(samplesKey, prev);
      },
    };
  },
  request: (raw) => {
    const p = flat<RemoveSampleTagInput, RemoveSampleTagScope>(raw);
    return api.removeSampleTag(p.sampleId, p.tagId, buildAuthOpts(p));
  },
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// add_tag (exposure)
// ---------------------------------------------------------------------------

export type AddExposureTagInput = { key: string; value: string };
type AddExposureTagScope = { sampleId: number; exposureId: number };

export const addExposureTagMutator: Mutator<OpPayload<AddExposureTagInput>, ExposureTag> = {
  kind: "add_tag",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<AddExposureTagInput, AddExposureTagScope>(raw);
    const placeholderId = nextOptimisticId();
    return {
      restore: rewriteExposureLists(qc, p.sampleId, (list) =>
        list.map((e) =>
          e.id === p.exposureId
            ? { ...e, tags: [...e.tags, {
                id: placeholderId, key: p.key, value: p.value, source: "manual",
              }] }
            : e,
        ),
      ),
    };
  },
  request: (raw) => {
    const p = flat<AddExposureTagInput, AddExposureTagScope>(raw);
    return api.addExposureTag(p.exposureId, p.key, p.value, buildAuthOpts(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flat<AddExposureTagInput, AddExposureTagScope>(raw);
    rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) => {
        if (e.id !== p.exposureId) return e;
        const filtered = e.tags.filter((t) =>
          !(t.id < 0 && t.key === p.key && t.value === p.value)
          && t.id !== response.id,
        );
        return { ...e, tags: [...filtered, response] };
      }),
    );
  },
};

// ---------------------------------------------------------------------------
// remove_tag (exposure)
// ---------------------------------------------------------------------------

export type RemoveExposureTagInput = { tagId: number };
type RemoveExposureTagScope = { sampleId: number; exposureId: number };

export const removeExposureTagMutator: Mutator<OpPayload<RemoveExposureTagInput>, void> = {
  kind: "remove_tag",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<RemoveExposureTagInput, RemoveExposureTagScope>(raw);
    return {
      restore: rewriteExposureLists(qc, p.sampleId, (list) =>
        list.map((e) =>
          e.id === p.exposureId
            ? { ...e, tags: e.tags.filter((t) => t.id !== p.tagId) }
            : e,
        ),
      ),
    };
  },
  request: (raw) => {
    const p = flat<RemoveExposureTagInput, RemoveExposureTagScope>(raw);
    return api.removeExposureTag(p.exposureId, p.tagId, buildAuthOpts(p));
  },
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// post_message
// ---------------------------------------------------------------------------

export type PostSampleMessageInput = { body: string };
type PostSampleMessageScope = { sampleId: number };

export const postSampleMessageMutator: Mutator<OpPayload<PostSampleMessageInput>, SampleMessage> = {
  kind: "post_message",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<PostSampleMessageInput, PostSampleMessageScope>(raw);
    const key = queryKeys.messages(p.sampleId);
    const prev = qc.getQueryData<SampleMessage[]>(key);
    const placeholder: SampleMessage = {
      id: nextOptimisticId(),
      sample_id: p.sampleId,
      author_id: null,
      author: p.username ?? null,
      body: p.body,
      created_at: new Date().toISOString(),
    };
    qc.setQueryData<SampleMessage[]>(key, [...(prev ?? []), placeholder]);
    return {
      restore: () => {
        if (prev === undefined) qc.removeQueries({ queryKey: key, exact: true });
        else qc.setQueryData(key, prev);
      },
    };
  },
  request: (raw) => {
    const p = flat<PostSampleMessageInput, PostSampleMessageScope>(raw);
    return api.postSampleMessage(p.sampleId, p.body, buildAuthOpts(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flat<PostSampleMessageInput, PostSampleMessageScope>(raw);
    const key = queryKeys.messages(p.sampleId);
    const list = qc.getQueryData<SampleMessage[]>(key) ?? [];
    // Replace the most recent negative-id placeholder for this body, and
    // dedupe against any concurrent SSE that already inserted the real msg.
    const seen = new Set<number>();
    const next: SampleMessage[] = [];
    let replaced = false;
    for (const m of list) {
      if (m.id < 0 && !replaced && m.body === response.body
          && m.sample_id === response.sample_id) {
        if (!seen.has(response.id)) { next.push(response); seen.add(response.id); }
        replaced = true;
        continue;
      }
      if (seen.has(m.id)) continue;
      next.push(m);
      seen.add(m.id);
    }
    if (!replaced && !seen.has(response.id)) next.push(response);
    qc.setQueryData<SampleMessage[]>(key, next);
  },
};

// ---------------------------------------------------------------------------
// set_exposure_status
// ---------------------------------------------------------------------------

export type SetExposureStatusInput = {
  exposureId: number;
  status: "accepted" | "rejected" | null;
};
type SetExposureStatusScope = { sampleId: number };

export const setExposureStatusMutator: Mutator<
  OpPayload<SetExposureStatusInput>, { id: number; status: string | null }
> = {
  kind: "set_exposure_status",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<SetExposureStatusInput, SetExposureStatusScope>(raw);
    const restoreExposures = rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) =>
        e.id === p.exposureId ? { ...e, status: p.status } : e,
      ),
    );
    const exposureKey = queryKeys.exposure(p.exposureId);
    const prevSingle = qc.getQueryData<Exposure>(exposureKey);
    if (prevSingle && prevSingle.id === p.exposureId) {
      qc.setQueryData<Exposure>(exposureKey, { ...prevSingle, status: p.status });
    }
    return {
      restore: () => {
        restoreExposures();
        if (prevSingle !== undefined) qc.setQueryData(exposureKey, prevSingle);
      },
    };
  },
  request: (raw) => {
    const p = flat<SetExposureStatusInput, SetExposureStatusScope>(raw);
    return api.setExposureStatus(p.exposureId, p.status, buildAuthOpts(p));
  },
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// select_exposure
// ---------------------------------------------------------------------------

export type SelectExposureInput = { exposureId: number };
type SelectExposureScope = { sampleId: number };

export const selectExposureMutator: Mutator<
  OpPayload<SelectExposureInput>, { id: number; selected: boolean }
> = {
  kind: "select_exposure",
  onMutate: (raw, qc): RollbackContext => {
    const p = flat<SelectExposureInput, SelectExposureScope>(raw);
    return {
      restore: rewriteExposureLists(qc, p.sampleId, (list) =>
        list.map((e) => ({ ...e, selected: e.id === p.exposureId })),
      ),
    };
  },
  request: (raw) => {
    const p = flat<SelectExposureInput, SelectExposureScope>(raw);
    return api.selectExposure(p.exposureId, buildAuthOpts(p));
  },
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Walk every cached `["sample", sampleId, "exposures", …]` query — multiple
 * variants exist (excludeRejected on/off) — apply `mutate` to each, and
 * return a `restore` closure that captures the snapshots for rollback.
 */
function rewriteExposureLists(
  qc: QueryClient,
  sampleId: number,
  mutate: (list: Exposure[]) => Exposure[],
): () => void {
  const prefix = ["sample", sampleId, "exposures"] as const;
  const entries = qc.getQueriesData<Exposure[]>({ queryKey: prefix });
  for (const [key, value] of entries) {
    if (!value) continue;
    qc.setQueryData<Exposure[]>(key, mutate(value));
  }
  return () => {
    for (const [key, value] of entries) {
      qc.setQueryData(key, value);
    }
  };
}
