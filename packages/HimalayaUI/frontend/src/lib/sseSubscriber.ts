import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../queries";

interface CurationEvent {
  kind?: string;
  entity_type?: string;
  entity_id?: number;
  actor?: string | null;
  client_id?: string | null;
}

/**
 * Handle a raw SSE "curation" event data string. Parses the JSON, applies the
 * self-echo filter, and invalidates TanStack Query caches for the affected
 * exposure.
 *
 * Extracted from App.tsx for unit-testability — callers pass a minimal ctx
 * object instead of the full React component context.
 */
export function handleCurationEvent(
  data: string,
  ctx: { clientId: string; qc: QueryClient },
): void {
  let event: CurationEvent | null = null;
  try {
    event = JSON.parse(data) as CurationEvent;
  } catch {
    return; // malformed frame, ignore
  }
  if (!event || typeof event.entity_id !== "number") return;
  // Self-echo filter: skip events authored by this tab. Other tabs of the
  // same user (or other users) pass through. System events (client_id=null)
  // also pass through — null !== string, so the equality check handles them
  // without a separate truthiness guard.
  //
  // Defense-in-depth: the queue framework's handleRemoteEvent (M1.2) covers
  // same-tab confirmation precisely via client_op_id matching for queue
  // mutations. This client_id-based filter is the broader fallback for the
  // legacy invalidation path that is still wired in App.tsx — it prevents
  // self-echo refetch for any event this tab originated, queue-shaped or
  // not. Safe to drop once App.tsx switches to handleRemoteEvent for all
  // SSE frames.
  if (event.client_id === ctx.clientId) return;
  if (event.entity_type !== "exposure") return;
  const id = event.entity_id;
  ctx.qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.exposure(id) });
}
