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
  Sample, Exposure,
  AuthOpts,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
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

export type UpdateSampleInput = { name?: string; notes?: string };
type UpdateSampleScope = BaseScope & { experimentId: number; sampleId: number };

export const updateSampleMutator: Mutator<UpdateSampleInput, UpdateSampleScope, Sample> = {
  kind: "update_sample",
  onMutate: (p, qc): RollbackContext => {
    const samplesKey = queryKeys.samples(p.experimentId);
    const prevList = qc.getQueryData<Sample[]>(samplesKey);
    const patch = patchOf(p);

    if (prevList) {
      qc.setQueryData<Sample[]>(samplesKey, prevList.map((s) =>
        s.id === p.sampleId ? { ...s, ...patch } : s,
      ));
    }
    return {
      restore: () => {
        if (prevList !== undefined) qc.setQueryData(samplesKey, prevList);
      },
    };
  },
  request: (p) => api.updateSample(p.sampleId, patchOf(p), buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // When SSE wins the race against HTTP, `synthesizeResponseFromSse`
    // hands us the SSE payload (the diff — only the patched field), so
    // unpatched fields are `undefined` on the response. Skip undefined
    // fields below so they don't clobber existing values.
    const samplesKey = queryKeys.samples(p.experimentId);
    const patch: Partial<Sample> = {};
    if (response.name !== undefined) patch.name = response.name;
    if (response.notes !== undefined) patch.notes = response.notes;
    const list = qc.getQueryData<Sample[]>(samplesKey);
    if (list) {
      qc.setQueryData<Sample[]>(samplesKey, list.map((s) =>
        s.id === p.sampleId ? { ...s, ...patch } : s));
    }
  },
};

function patchOf(p: { name?: string; notes?: string }): { name?: string; notes?: string } {
  const out: { name?: string; notes?: string } = {};
  if (p.name !== undefined) out.name = p.name;
  if (p.notes !== undefined) out.notes = p.notes;
  return out;
}

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
    return { restore: restoreExposures };
  },
  request: (p) => api.setExposureStatus(p.exposureId, p.status, buildAuthOpts(p)),
  onSuccess: (_p, _r, qc) => {
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
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
  onSuccess: (_p, _r, qc) => {
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
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
