import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import type { Trace } from "../api";
import { queryKeys } from "../queries";
import { representativeExposure } from "../lib/exposure";

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
  });
  const out = new Map<number, Trace>();
  sampleIds.forEach((id, i) => {
    const t = traces[i]?.data;
    if (t) out.set(id, t);
  });
  return out;
}
