import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../queries";

interface CurationEvent {
  kind?: string;
  entity_type?: string;
  entity_id?: number;
  actor?: string | null;
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
  ctx: { username: string | undefined; qc: QueryClient },
): void {
  let event: CurationEvent | null = null;
  try {
    event = JSON.parse(data) as CurationEvent;
  } catch {
    return; // malformed frame, ignore
  }
  if (!event || typeof event.entity_id !== "number") return;
  // Self-echo filter: skip events authored by this client. If two tabs share
  // the same X-Username (same lab user), both tabs' echoes are filtered.
  if (event.actor && event.actor === ctx.username) return;
  if (event.entity_type !== "exposure") return;
  const id = event.entity_id;
  ctx.qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
  ctx.qc.invalidateQueries({ queryKey: queryKeys.exposure(id) });
}
