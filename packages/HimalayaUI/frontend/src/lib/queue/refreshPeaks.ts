import type { Mutation, QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../queries";

const PEAK_KINDS = new Set(["peak_added", "peak_excluded", "peak_unexcluded", "peak_removed"]);

const isPendingPeakOp = (m: Mutation, exposureId: number) =>
  m.state.status === "pending" &&
  PEAK_KINDS.has(String(m.options.mutationKey?.[0])) &&
  (m.state.variables as { exposureId?: number } | undefined)?.exposureId === exposureId;

const waiting = new WeakMap<QueryClient, Set<number>>();

export function refreshPeaks(qc: QueryClient, exposureId: number): void {
  const mc = qc.getMutationCache();
  const run = () => void qc.invalidateQueries({ queryKey: queryKeys.peaks(exposureId) });
  const blocked = () => mc.getAll().some((m) => isPendingPeakOp(m, exposureId));
  if (!blocked()) return run();
  const set = waiting.get(qc) ?? new Set<number>();
  waiting.set(qc, set);
  if (set.has(exposureId)) return;
  set.add(exposureId);
  const unsubscribe = mc.subscribe(() => {
    if (blocked()) return;
    unsubscribe();
    set.delete(exposureId);
    run();
  });
}

export function cancelPeakRefetch(qc: QueryClient, exposureId: number): void {
  void qc.cancelQueries({ queryKey: queryKeys.peaks(exposureId) });
}
