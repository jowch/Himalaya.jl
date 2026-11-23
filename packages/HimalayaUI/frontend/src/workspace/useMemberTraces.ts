import { useMemo } from "react";
import { useQueries } from "@tanstack/react-query";
import * as api from "../api";
import type { Trace } from "../api";
import { queryKeys, retryUnlessMissing } from "../queries";
import type { MemberLine } from "../ui/WorkspacePlot";

function dataOf(results: readonly { data?: Trace | undefined }[]): (Trace | undefined)[] {
  return results.map((r) => r.data);
}

export function useMemberTraces(exposureIds: readonly number[]): MemberLine[] {
  const traces = useQueries({
    queries: exposureIds.map((id) => ({
      queryKey: queryKeys.trace(id),
      queryFn: () => api.getTrace(id),
      retry: retryUnlessMissing,
    })),
    combine: dataOf,
  });
  const key = exposureIds.join(",");
  return useMemo(() => {
    const ids = key ? key.split(",").map(Number) : [];
    return ids.flatMap((id, i) => {
      const t = traces[i];
      return t && t.q.length > 1 ? [{ id, q: t.q, I: t.I }] : [];
    });
  }, [key, traces]);
}
