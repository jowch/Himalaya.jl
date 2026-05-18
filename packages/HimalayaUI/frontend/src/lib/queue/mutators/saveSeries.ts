/**
 * saveSeries mutator (#166). One mutator handles BOTH create and recipe edit;
 * `payload.id` discriminates — absent ⇒ create (POST /api/series), present ⇒
 * recipe edit (PATCH /api/series/:id). Mirrors `saveComparison`.
 *
 * No optimistic write — the builder UI (I3.5b) shows the local draft; this
 * mutator's job is the server round-trip + cache reconciliation. `onSuccess`
 * splices a full Series response into the detail cache and invalidates the
 * listing; on the SSE-wins path the response is a partial shape, detected and
 * routed to invalidate-fallback.
 */
import * as api from "../../../api";
import type { AuthOpts, Series, SaveSeriesBody, SeriesSampleInput } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface SaveSeriesInput {
  /** Absent ⇒ create; present ⇒ recipe edit of an existing series. */
  id?: number;
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  ordering_variable?: string | null;
  order_rule?: "ascending" | "descending" | "manual";
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
  view_grouping_mode?: string | null;
  view_show_peak_ticks?: boolean | null;
  view_show_peak_labels?: boolean | null;
}

interface SaveSeriesScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

function buildBody(p: SaveSeriesInput): SaveSeriesBody {
  const out: SaveSeriesBody = { title: p.title, samples: p.samples };
  if (p.description !== undefined) out.description = p.description;
  if (p.ordering_variable !== undefined) out.ordering_variable = p.ordering_variable;
  if (p.order_rule !== undefined) out.order_rule = p.order_rule;
  if (p.forked_from_id !== undefined) out.forked_from_id = p.forked_from_id;
  if (p.forked_at_hash !== undefined) out.forked_at_hash = p.forked_at_hash;
  if (p.view_grouping_mode !== undefined) out.view_grouping_mode = p.view_grouping_mode;
  if (p.view_show_peak_ticks !== undefined) out.view_show_peak_ticks = p.view_show_peak_ticks;
  if (p.view_show_peak_labels !== undefined) out.view_show_peak_labels = p.view_show_peak_labels;
  return out;
}

export const saveSeriesMutator: Mutator<SaveSeriesInput, SaveSeriesScope, Series> = {
  kind: "series_save",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.saveSeries(buildBody(p), p.id, buildAuthOpts(p)),
  onSuccess: (_p, response, qc) => {
    // SSE-wins path: `synthesizeFromSse` yields a partial shape (no `samples`
    // array, no `state`). Probe for the full Series shape; fall back to
    // invalidate so the next read fetches the canonical projection.
    const looksFull = Array.isArray((response as { samples?: unknown }).samples)
      && typeof (response as { state?: unknown }).state === "string";
    if (looksFull) {
      qc.setQueryData(queryKeys.series(response.id), response);
    } else if (typeof response?.id === "number") {
      qc.invalidateQueries({ queryKey: queryKeys.series(response.id) });
    }
    qc.invalidateQueries({ queryKey: queryKeys.seriesList });
  },
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    // Partial Series shape — `onSuccess`'s looksFull detector trips the
    // invalidate fallback (no `state` in the SSE payload). `id` is required so
    // the invalidate branch targets the right cache key.
    return { ...base, ...payload, id: remote.entity_id } as unknown as Series;
  },
};
