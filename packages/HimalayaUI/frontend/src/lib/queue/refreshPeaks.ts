import type { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../queries";

export function refreshPeaks(qc: QueryClient, exposureId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.peaks(exposureId) });
}
