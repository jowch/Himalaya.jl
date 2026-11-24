import { useCallback, useEffect, useMemo, useState } from "react";
import { useQueryClient } from "@tanstack/react-query";
import * as api from "../api";
import { queryKeys } from "../queries";

const DESCRIPTION_FRESH_MS = 30_000;

export interface FolderDescriptions {
  found: ReadonlyMap<string, api.FolderDescription>;
  describeNow: (path: string) => Promise<api.FolderDescription>;
}

export function useFolderDescriptions(paths: ReadonlyArray<string>, enabled: boolean): FolderDescriptions {
  const qc = useQueryClient();
  const joined = useMemo(() => paths.join("\n"), [paths]);
  const [found, setFound] = useState<ReadonlyMap<string, api.FolderDescription>>(new Map());

  const describe = useCallback(
    async (path: string, signal?: AbortSignal): Promise<api.FolderDescription> => {
      const d = await api.describeFolder(path, signal === undefined ? undefined : { signal });
      qc.setQueryData(queryKeys.folderDescription(path), d);
      return d;
    },
    [qc],
  );

  useEffect(() => {
    const wanted = joined === "" ? [] : joined.split("\n");
    if (!enabled) { setFound(new Map()); return; }
    const cached = new Map<string, api.FolderDescription>();
    for (const p of wanted) {
      const state = qc.getQueryState<api.FolderDescription>(queryKeys.folderDescription(p));
      if (state?.data !== undefined && !state.isInvalidated && Date.now() - state.dataUpdatedAt < DESCRIPTION_FRESH_MS) cached.set(p, state.data);
    }
    setFound(cached);
    const ctrl = new AbortController();
    void (async () => {
      for (const p of wanted) {
        if (ctrl.signal.aborted) return;
        if (cached.has(p)) continue;
        try {
          const d = await describe(p, ctrl.signal);
          if (ctrl.signal.aborted) return;
          setFound((m) => new Map(m).set(p, d));
        } catch {
          if (ctrl.signal.aborted) return;
        }
      }
    })();
    return () => ctrl.abort();
  }, [joined, enabled, qc, describe]);

  const describeNow = useCallback((path: string) => describe(path), [describe]);
  return { found, describeNow };
}
