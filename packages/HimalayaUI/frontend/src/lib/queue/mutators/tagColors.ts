import * as api from "../../../api";
import type { TagColor, TagColorEntry } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { applyTagColor } from "../../labels";
import type { Mutator } from "../types";

interface TagColorScope {
  username: string | undefined;
  clientId: string;
}

type TagColorInput = { label: string; color: TagColor | null };

export const setTagColorMutator: Mutator<TagColorInput, TagColorScope, TagColorEntry[]> = {
  kind: "set_tag_color",
  onMutate: (p, qc) => {
    const key = queryKeys.tagColors();
    const prev = qc.getQueryData<TagColorEntry[]>(key);
    qc.setQueryData<TagColorEntry[]>(key, applyTagColor(prev ?? [], p.label, p.color));
    return { restore: () => qc.setQueryData(key, prev) };
  },
  request: (p) => api.setTagColor(p.label, p.color, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (_p, response, qc) => {
    if (Array.isArray(response)) qc.setQueryData<TagColorEntry[]>(queryKeys.tagColors(), response as TagColorEntry[]);
    qc.invalidateQueries({ queryKey: queryKeys.tagColors() });
  },
};
