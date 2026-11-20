import type { MutationCache, Query, QueryClient } from "@tanstack/react-query";

const isRawTrace = (q: Query) => q.queryKey[0] === "exposure" && q.queryKey[2] === "trace";

const hasPending = (mc: MutationCache) => mc.getAll().some((m) => m.state.status === "pending");

export function createReconnectRefresh(qc: QueryClient, mc: MutationCache): () => void {
  let scheduled = false;
  const run = () => {
    scheduled = false;
    void qc.invalidateQueries({ predicate: (q) => !isRawTrace(q) });
  };
  return () => {
    if (scheduled) return;
    scheduled = true;
    if (!hasPending(mc)) return run();
    const unsubscribe = mc.subscribe(() => {
      if (hasPending(mc)) return;
      unsubscribe();
      run();
    });
  };
}
