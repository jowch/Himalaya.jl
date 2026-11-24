import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import type { Trace } from "../api";
import { queryKeys, retryUnlessMissing } from "../queries";
import { analyzedExposure } from "../lib/exposure";

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
  const reps = exposures.map((q) => analyzedExposure(q.data ?? [])?.id);
  const traces = useQueries({
    queries: reps.map((eid) => ({
      queryKey: queryKeys.trace(eid),
      queryFn: () => api.getTrace(eid as number),
      enabled: eid !== undefined,
      retry: retryUnlessMissing,
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

export type TraceStatus = "loading" | "missing";

export function useTraceStatuses(sampleIds: readonly number[]): Map<number, TraceStatus> {
  const exposures = useQueries({
    queries: sampleIds.map((id) => ({
      queryKey: queryKeys.exposures(id),
      queryFn: () => api.listExposures(id),
    })),
  });
  const reps = exposures.map((q) => analyzedExposure(q.data ?? [])?.id);
  const fetching = useQueries({
    queries: reps.map((eid) => ({
      queryKey: queryKeys.trace(eid),
      queryFn: () => api.getTrace(eid as number),
      enabled: eid !== undefined,
      retry: retryUnlessMissing,
    })),
    combine: (results) => results.map((r) => r.isPending && r.fetchStatus !== "idle"),
  });
  const key = sampleIds.join(",");
  const flags = exposures.map((q, i) => q.isPending || fetching[i] === true).join(",");
  return useMemo(() => {
    const out = new Map<number, TraceStatus>();
    const loading = flags.split(",");
    key.split(",").forEach((id, i) => {
      if (id) out.set(Number(id), loading[i] === "true" ? "loading" : "missing");
    });
    return out;
  }, [key, flags]);
}
