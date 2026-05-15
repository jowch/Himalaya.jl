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
  Sample, SampleTag, Exposure, ExposureTag, SampleMessage, ComparisonMessage,
  AuthOpts,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { nextOptimisticId } from "../optimisticId";
import { replacePlaceholder } from "../replacePlaceholder";
import type { Mutator, RollbackContext } from "../types";

interface BaseScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: {
  username?: string | undefined;
  clientId?: string | undefined;
  clientOpId?: string | undefined;
}): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

// ---------------------------------------------------------------------------
// update_sample
// ---------------------------------------------------------------------------

export type UpdateSampleInput = { display_name?: string; notes?: string };
type UpdateSampleScope = BaseScope & { experimentId: number; sampleId: number };

export const updateSampleMutator: Mutator<UpdateSampleInput, UpdateSampleScope, Sample> = {
  kind: "update_sample",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.updateSample(p.sampleId, patchOf(p), buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // The PATCH route returns only the `samples` row — no `tags` field.
    // The Sample type requires `tags: SampleTag[]`, and the listing routes
    // (GET /api/experiments/:id/samples and GET /api/samples/:id) attach
    // tags client-side. Spreading the response wholesale into the cache
    // would clobber `tags` to undefined.  Merge ONLY the patched fields
    // onto the existing cache entry so tags survive.
    //
    // Also: when SSE wins the race against HTTP, `synthesizeResponseFromSse`
    // hands us the SSE payload (the diff — only the patched field), so
    // unpatched fields are `undefined` on the response. Skip undefined
    // fields below so they don't clobber existing values.
    const samplesKey = queryKeys.samples(p.experimentId);
    const sampleKey = queryKeys.sample(p.sampleId);
    const patch: Partial<Sample> = {};
    if (response.name !== undefined) patch.name = response.name;
    if (response.notes !== undefined) patch.notes = response.notes;
    if (response.display_name !== undefined) patch.display_name = response.display_name;
    const list = qc.getQueryData<Sample[]>(samplesKey);
    if (list) {
      qc.setQueryData<Sample[]>(samplesKey, list.map((s) =>
        s.id === p.sampleId ? { ...s, ...patch } : s));
    }
    const prevSingle = qc.getQueryData<Sample>(sampleKey);
    if (prevSingle && prevSingle.id === p.sampleId) {
      qc.setQueryData<Sample>(sampleKey, { ...prevSingle, ...patch });
    }
  },
};

function patchOf(p: { display_name?: string; notes?: string }): { display_name?: string; notes?: string } {
  const out: { display_name?: string; notes?: string } = {};
  if (p.display_name !== undefined) out.display_name = p.display_name;
  if (p.notes !== undefined) out.notes = p.notes;
  return out;
}

// ---------------------------------------------------------------------------
// add_tag (sample)
// ---------------------------------------------------------------------------

export type AddSampleTagInput = { key: string; value: string };
type AddSampleTagScope = BaseScope & { experimentId: number; sampleId: number };

export const addSampleTagMutator: Mutator<AddSampleTagInput, AddSampleTagScope, SampleTag> = {
  kind: "add_tag",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.addSampleTag(p.sampleId, p.key, p.value, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // The route emits `{id, sample_id, key, value, source}` (routes_samples.jl)
    // but SampleTag omits `sample_id`. Strip it so the cache row matches type.
    const tag: SampleTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    const samplesKey = queryKeys.samples(p.experimentId);
    qc.setQueryData<Sample[]>(samplesKey, (list) => {
      if (!list) return list;
      return list.map((s) =>
        s.id !== p.sampleId
          ? s
          : {
              ...s,
              tags: replacePlaceholder(
                s.tags,
                tag,
                (t) => t.key === p.key && t.value === p.value,
              ),
            },
      );
    });
  },
};

// ---------------------------------------------------------------------------
// remove_tag (sample)
// ---------------------------------------------------------------------------

export type RemoveSampleTagInput = { tagId: number };
type RemoveSampleTagScope = BaseScope & { experimentId: number; sampleId: number };

export const removeSampleTagMutator: Mutator<RemoveSampleTagInput, RemoveSampleTagScope, void> = {
  kind: "remove_tag",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.removeSampleTag(p.sampleId, p.tagId, buildAuthOpts(p)),
  onSuccess: () => {},
  // 404 = the tag is already gone → desired end state.
  treats404AsSuccess: true,
};

// ---------------------------------------------------------------------------
// add_tag (exposure)
// ---------------------------------------------------------------------------

export type AddExposureTagInput = { key: string; value: string };
type AddExposureTagScope = BaseScope & { sampleId: number; exposureId: number };

export const addExposureTagMutator: Mutator<AddExposureTagInput, AddExposureTagScope, ExposureTag> = {
  kind: "add_tag",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.addExposureTag(p.exposureId, p.key, p.value, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // The route emits `{id, exposure_id, key, value, source}` (routes_exposures.jl)
    // but ExposureTag omits `exposure_id`. Strip it before caching.
    const tag: ExposureTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) =>
        e.id !== p.exposureId
          ? e
          : {
              ...e,
              tags: replacePlaceholder(
                e.tags,
                tag,
                (t) => t.key === p.key && t.value === p.value,
              ),
            },
      ),
    );
  },
};

// ---------------------------------------------------------------------------
// remove_tag (exposure)
// ---------------------------------------------------------------------------

export type RemoveExposureTagInput = { tagId: number };
type RemoveExposureTagScope = BaseScope & { sampleId: number; exposureId: number };

export const removeExposureTagMutator: Mutator<RemoveExposureTagInput, RemoveExposureTagScope, void> = {
  kind: "remove_tag",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.removeExposureTag(p.exposureId, p.tagId, buildAuthOpts(p)),
  onSuccess: () => {},
  // 404 = the tag is already gone → desired end state.
  treats404AsSuccess: true,
};

// ---------------------------------------------------------------------------
// post_message
// ---------------------------------------------------------------------------

export type PostSampleMessageInput = { body: string };
type PostSampleMessageScope = BaseScope & { sampleId: number };

export const postSampleMessageMutator: Mutator<PostSampleMessageInput, PostSampleMessageScope, SampleMessage> = {
  kind: "post_message",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.postSampleMessage(p.sampleId, p.body, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    const key = queryKeys.messages(p.sampleId);
    qc.setQueryData<SampleMessage[]>(key, (list) =>
      replacePlaceholder(
        list ?? [],
        response,
        (m) => m.body === response.body && m.sample_id === response.sample_id,
      ),
    );
  },
};

// ---------------------------------------------------------------------------
// post_message — comparison context (Phase 10)
// ---------------------------------------------------------------------------
// Same OpKind ("post_message") as the sample variant; the registry
// discriminates by payload shape (`comparisonId` vs `sampleId`). Cache key
// is `queryKeys.comparisonMessages(comparisonId)` and the wire route is
// `POST /api/comparisons/:id/messages`. The placeholder shape mirrors the
// `ComparisonMessage` row.

export type PostComparisonMessageInput = { body: string };
type PostComparisonMessageScope = BaseScope & { comparisonId: number };

export const postComparisonMessageMutator: Mutator<
  PostComparisonMessageInput, PostComparisonMessageScope, ComparisonMessage
> = {
  kind: "post_message",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.comparisonMessages(p.comparisonId);
    const prev = qc.getQueryData<ComparisonMessage[]>(key);
    const placeholder: ComparisonMessage = {
      id: nextOptimisticId(),
      comparison_id: p.comparisonId,
      author_id: null,
      author: p.username ?? null,
      body: p.body,
      created_at: new Date().toISOString(),
    };
    qc.setQueryData<ComparisonMessage[]>(key, [...(prev ?? []), placeholder]);
    return {
      restore: () => {
        if (prev === undefined) qc.removeQueries({ queryKey: key, exact: true });
        else qc.setQueryData(key, prev);
      },
    };
  },
  request: (p) => api.postComparisonMessage(p.comparisonId, p.body, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    const key = queryKeys.comparisonMessages(p.comparisonId);
    qc.setQueryData<ComparisonMessage[]>(key, (list) =>
      replacePlaceholder(
        list ?? [],
        response,
        (m) => m.body === response.body && m.comparison_id === response.comparison_id,
      ),
    );
  },
};

// ---------------------------------------------------------------------------
// set_exposure_status
// ---------------------------------------------------------------------------

export type SetExposureStatusInput = {
  exposureId: number;
  status: "accepted" | "rejected" | null;
};
type SetExposureStatusScope = BaseScope & { sampleId: number };

export const setExposureStatusMutator: Mutator<
  SetExposureStatusInput, SetExposureStatusScope, { id: number; status: string | null }
> = {
  kind: "set_exposure_status",
  onMutate: (p, qc): RollbackContext => {
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
  request: (p) => api.setExposureStatus(p.exposureId, p.status, buildAuthOpts(p)),
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// select_exposure
// ---------------------------------------------------------------------------

export type SelectExposureInput = { exposureId: number };
type SelectExposureScope = BaseScope & { sampleId: number };

export const selectExposureMutator: Mutator<
  SelectExposureInput, SelectExposureScope, { id: number; selected: boolean }
> = {
  kind: "select_exposure",
  onMutate: (p, qc): RollbackContext => {
    return {
      restore: rewriteExposureLists(qc, p.sampleId, (list) =>
        list.map((e) => ({ ...e, selected: e.id === p.exposureId })),
      ),
    };
  },
  request: (p) => api.selectExposure(p.exposureId, buildAuthOpts(p)),
  onSuccess: () => {},
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Walk every cached `["sample", sampleId, "exposures", …]` query, apply
 * `mutate` to each, and return a `restore` closure that captures the
 * snapshots for rollback. Prefix-walks rather than direct-keys so any
 * future per-view variant still gets the rewrite without extra plumbing.
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
