import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { Experiment, Setup } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator } from "../types";

interface SetupScope {
  experimentId: number;
  username: string | undefined;
  clientId: string;
}

export type SetupPatch = { name?: string | null; included?: boolean };

export type SetupFields = Partial<Pick<Setup, "name" | "display_name" | "included">>;

function withSetup(e: Experiment, setupId: number, patch: SetupFields): Experiment {
  if (!e.setups) return e;
  return { ...e, setups: e.setups.map((s) => (s.id === setupId ? { ...s, ...patch } : s)) };
}

export function mergeSetupPatch(qc: QueryClient, experimentId: number, setupId: number, patch: SetupFields): void {
  qc.setQueryData<Experiment>(queryKeys.experiment(experimentId), (cur) => (cur ? withSetup(cur, setupId, patch) : cur));
  qc.setQueryData<Experiment[]>(queryKeys.experiments, (list) =>
    list?.map((e) => (e.id === experimentId ? withSetup(e, setupId, patch) : e)));
}

function fieldsOf(qc: QueryClient, experimentId: number, setupId: number, p: SetupPatch): SetupFields {
  const out: SetupFields = {};
  if (p.included !== undefined) out.included = p.included;
  if (p.name !== undefined) {
    const folder = qc.getQueryData<Experiment>(queryKeys.experiment(experimentId))?.setups?.find((s) => s.id === setupId)?.folder;
    const typed = p.name === null ? "" : p.name.trim();
    out.display_name = typed === "" || typed === folder ? null : typed;
    out.name = out.display_name ?? folder ?? typed;
  }
  return out;
}

function requestOf(p: SetupPatch): SetupPatch {
  const out: SetupPatch = {};
  if (p.name !== undefined) out.name = p.name;
  if (p.included !== undefined) out.included = p.included;
  return out;
}

export const patchSetupMutator: Mutator<SetupPatch & { setupId: number }, SetupScope, Setup> = {
  kind: "patch_setup",
  onMutate: (p, qc) => {
    const prevOne = qc.getQueryData<Experiment>(queryKeys.experiment(p.experimentId));
    const prevList = qc.getQueryData<Experiment[]>(queryKeys.experiments);
    mergeSetupPatch(qc, p.experimentId, p.setupId, fieldsOf(qc, p.experimentId, p.setupId, p));
    return {
      restore: () => {
        qc.setQueryData(queryKeys.experiment(p.experimentId), prevOne);
        qc.setQueryData(queryKeys.experiments, prevList);
      },
    };
  },
  request: (p) => api.patchSetup(p.setupId, requestOf(p), authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => {
    const patch: SetupFields = {};
    if (response.name !== undefined) patch.name = response.name;
    if (response.display_name !== undefined) patch.display_name = response.display_name;
    if (response.included !== undefined) patch.included = response.included;
    mergeSetupPatch(qc, p.experimentId, p.setupId, patch);
  },
};
