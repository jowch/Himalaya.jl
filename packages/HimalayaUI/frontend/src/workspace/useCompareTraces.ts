import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import type { Trace } from "../api";
import { queryKeys } from "../queries";
import { representativeExposure } from "../lib/exposure";

function dataOf(results: readonly { data?: Trace | undefined }[]): (Trace | undefined)[] {
  return results.map((r) => r.data);
}

export function useCompareTraces(sampleIds: readonly number[]): Map<number, Trace> {
  const exposures = useQueries({
    queries: sampleIds.map((id) => ({
      queryKey: queryKeys.exposures(id),
      queryFn: () => api.listExposures(id),
    })),
  });
  const reps = exposures.map((q) => representativeExposure(q.data ?? [])?.id);
  const traces = useQueries({
    queries: reps.map((eid) => ({
      queryKey: queryKeys.trace(eid),
      queryFn: () => api.getTrace(eid as number),
      enabled: eid !== undefined,
    })),
    combine: dataOf,
  });
  const key = sampleIds.join(",");
  return useMemo(() => {
    const out = new Map<number, Trace>();
    key.split(",").forEach((id, i) => {
      const t = traces[i];
      if (id && t) out.set(Number(id), t);
    });
    return out;
  }, [key, traces]);
}
