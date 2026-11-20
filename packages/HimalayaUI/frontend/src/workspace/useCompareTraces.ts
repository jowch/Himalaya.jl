import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import type { Exposure, Trace } from "../api";
import { queryKeys } from "../queries";
import { representativeExposure } from "../lib/exposure";

export interface CompareRef {
  trace: Trace;
  exposure: Exposure;
}

function dataOf(results: readonly { data?: Trace | undefined }[]): (Trace | undefined)[] {
  return results.map((r) => r.data);
}

function representatives(results: readonly { data?: Exposure[] | undefined }[]): (Exposure | undefined)[] {
  return results.map((r) => representativeExposure(r.data ?? []));
}

export function useCompareTraces(sampleIds: readonly number[]): Map<number, CompareRef> {
  const reps = useQueries({
    queries: sampleIds.map((id) => ({
      queryKey: queryKeys.exposures(id),
      queryFn: () => api.listExposures(id),
    })),
    combine: representatives,
  });
  const traces = useQueries({
    queries: reps.map((e) => ({
      queryKey: queryKeys.trace(e?.id),
      queryFn: () => api.getTrace(e?.id as number),
      enabled: e !== undefined,
    })),
    combine: dataOf,
  });
  const key = sampleIds.join(",");
  return useMemo(() => {
    const out = new Map<number, CompareRef>();
    key.split(",").forEach((id, i) => {
      const t = traces[i];
      const e = reps[i];
      if (id && t && e) out.set(Number(id), { trace: t, exposure: e });
    });
    return out;
  }, [key, traces, reps]);
}
