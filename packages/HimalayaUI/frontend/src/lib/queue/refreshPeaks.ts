import type { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../queries";

/**
 * Refetch an exposure's peak list. The server decides which auto peaks an
 * exclude hides or an add replaces (grid-aware, nearest-first), so cache
 * patches stay for responsiveness and this refetch settles them.
 */
export function refreshPeaks(qc: QueryClient, exposureId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.peaks(exposureId) });
}
