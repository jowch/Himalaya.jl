import { refetchListsWith, refetchSelection } from "../selectionLists";
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
export type SetExposureStatusInput = {
  exposureId: number;
  status: "accepted" | "rejected" | null;
  selected?: boolean;
};
type SetExposureStatusScope = BaseScope & { sampleId: number };

export const setExposureStatusMutator: Mutator<
  SetExposureStatusInput, SetExposureStatusScope, { id: number; status: string | null }
> = {
  kind: "set_exposure_status",
  onMutate: (p, qc): RollbackContext => {
    const restoreExposures = rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) =>
        e.id === p.exposureId ? { ...e, status: p.status, ...(p.selected === undefined ? {} : { selected: p.selected }) } : e,
      ),
    );
    return { restore: restoreExposures };
  },
  request: (p) => api.setExposureStatus(p.exposureId, p.status, buildAuthOpts(p), p.selected),
  onSuccess: (p, _r, qc) => {
    refetchListsWith(qc, p.exposureId);
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
};
export type SetSelectionInput = { exposureIds: number[] };
type SetSelectionScope = BaseScope & { sampleId: number };

export const setSelectionMutator: Mutator<
  SetSelectionInput, SetSelectionScope,
  { sample_id: number; exposure_ids: number[]; analyzed_exposure_id: number | null }
> = {
  kind: "set_selection",
  onMutate: (p, qc): RollbackContext => {
    const chosen = new Set(p.exposureIds);
    const restoreExposures = rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) => (e.kind === "combined" || e.sample_id !== p.sampleId ? e : { ...e, selected: chosen.has(e.id) })),
    );
    return { restore: restoreExposures };
  },
  request: (p) => api.setSelection(p.sampleId, p.exposureIds, buildAuthOpts(p)),
  onSuccess: (p, _r, qc) => refetchSelection(qc, p.sampleId),
};
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