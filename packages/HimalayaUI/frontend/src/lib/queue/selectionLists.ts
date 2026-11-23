import type { QueryClient } from "@tanstack/react-query";
import type { Exposure } from "../../api";
import { isCombined } from "../exposure";

export function refetchListsWithoutCombined(qc: QueryClient, exposureId: number): void {
  const lists = qc.getQueriesData<Exposure[]>({
    predicate: (q) => q.queryKey[0] === "sample" && q.queryKey[2] === "exposures",
  });
  for (const [queryKey, list] of lists) {
    if (list?.some((e) => e.id === exposureId) && !list.some(isCombined)) {
      void qc.invalidateQueries({ queryKey, exact: true });
    }
  }
}
