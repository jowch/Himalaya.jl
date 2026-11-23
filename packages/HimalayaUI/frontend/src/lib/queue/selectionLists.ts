import type { Query, QueryClient } from "@tanstack/react-query";
import type { Exposure } from "../../api";
import { queryKeys } from "../../queries";

export function refetchListsWith(qc: QueryClient, exposureId: number): void {
  const lists = qc.getQueriesData<Exposure[]>({
    predicate: (q) => q.queryKey[0] === "sample" && q.queryKey[2] === "exposures",
  });
  for (const [queryKey, list] of lists) {
    if (list?.some((e) => e.id === exposureId)) void qc.invalidateQueries({ queryKey, exact: true });
  }
}

export function refetchSelection(qc: QueryClient, sampleId: number): void {
  void qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures"] });
  void qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  void qc.invalidateQueries({ predicate: (q) => q.queryKey[0] === "experiment" && q.queryKey[2] === "samples" });
}

const LIST_COALESCE_MS = 250;
const pendingListRefetch = new WeakMap<QueryClient, ReturnType<typeof setTimeout>>();

export function refetchAverage(qc: QueryClient, sampleId: number): void {
  void qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures"] });
  if (pendingListRefetch.has(qc)) return;
  pendingListRefetch.set(qc, setTimeout(() => {
    pendingListRefetch.delete(qc);
    void qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
    void qc.invalidateQueries({ predicate: (q) => q.queryKey[0] === "experiment" && q.queryKey[2] === "samples" });
  }, LIST_COALESCE_MS));
}

export function patchListsWith<T extends { id: number }>(
  qc: QueryClient,
  predicate: (q: Query) => boolean,
  id: number,
  patch: (row: T) => T,
): void {
  for (const [queryKey, list] of qc.getQueriesData<T[]>({ predicate })) {
    if (list?.some((row) => row.id === id)) qc.setQueryData<T[]>(queryKey, list.map((row) => (row.id === id ? patch(row) : row)));
  }
}
