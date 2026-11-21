import * as api from "../../../api";
import type { LabelColor, LabelColorEntry } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { applyLabelColor } from "../../labels";
import type { Mutator } from "../types";

interface LabelColorScope {
  username: string | undefined;
  clientId: string;
}

type LabelColorInput = { label: string; color: LabelColor | null };

export const setLabelColorMutator: Mutator<LabelColorInput, LabelColorScope, LabelColorEntry[]> = {
  kind: "set_tag_color",
  onMutate: (p, qc) => {
    const key = queryKeys.labelColors();
    const prev = qc.getQueryData<LabelColorEntry[]>(key);
    qc.setQueryData<LabelColorEntry[]>(key, applyLabelColor(prev ?? [], p.label, p.color));
    return { restore: () => qc.setQueryData(key, prev) };
  },
  request: (p) => api.setLabelColor(p.label, p.color, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (_p, response, qc) => {
    if (Array.isArray(response)) qc.setQueryData<LabelColorEntry[]>(queryKeys.labelColors(), response as LabelColorEntry[]);
    qc.invalidateQueries({ queryKey: queryKeys.labelColors() });
  },
};
