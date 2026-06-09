import type { Mutator, RollbackContext } from "../types";
import type { AuthOpts, Series, SeriesSampleInput } from "../../../api";
import { authOpts } from "../../authOpts";
import { queryKeys } from "../../../queries";
import * as api from "../../../api";

export interface ScopeAndCreateSeriesInput {
  /** The ordering-variable key to write to sample_tags (source='scoping'). */
  key: string;
  tags: { sampleId: number; value: string }[];
  /** SaveSeriesBody fields for the new series. */
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  orderingVariable?: string | null;
}

interface ScopeAndCreateSeriesScope {
  username: string | undefined;
  clientId: string;
}

// Local auth-opts wrapper — same pattern every mutator file uses (saveSeries,
// scopeSeries, …): wrap the shared authOpts(username, clientId, clientOpId) so
// request() reads the op's metadata off the flat payload.
function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/**
 * scopeAndCreateSeriesMutator — atomic scope-then-create (M-A Task 7).
 *
 * `request` runs sequentially:
 *   1. batchSampleTags (source='scoping') — writes the ordering tags FIRST.
 *   2. saveSeries (no id → POST /api/series) — then creates the series.
 * Returns the created `Series`.
 *
 * Uses kind='series_save' so the SSE frame the backend emits for the create
 * routes through the existing `saveSeriesMutator.synthesizeFromSse` path for
 * foreign-tab replay (mutatorRegistry maps the `series_save` event kind →
 * saveSeriesMutator). The scoping batch tags emit `add_tag` frames; those
 * replay via the existing `addSampleTagMutator.synthesizeFromSse` path. So
 * this mutator defines NO `synthesizeFromSse` of its own — it would be dead
 * code (the registry never resolves a foreign event to this mutator), exactly
 * as scopeSeriesMutator documents.
 *
 * onMutate is a no-op (mirrors scopeSeriesMutator): on success the user
 * navigates away to the new builder, so there is no on-surface cache to
 * optimistically patch/roll back.
 *
 * onSuccess: splices the new Series into the listing cache and primes the
 * detail cache (so the navigation to /series/:id reads it immediately), then
 * invalidates the two scoping proposal sources (same as scopeSeriesMutator).
 */
export const scopeAndCreateSeriesMutator: Mutator<
  ScopeAndCreateSeriesInput,
  ScopeAndCreateSeriesScope,
  Series
> = {
  kind: "series_save",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: async (p, _signal) => {
    const opts = buildAuthOpts(p);
    // 1. Write the ordering tags FIRST.
    await api.batchSampleTags(
      p.key,
      p.tags.map((t) => ({ sample_id: t.sampleId, value: t.value })),
      "scoping",
      opts,
    );
    // 2. Then create the series (no id → POST /api/series).
    return api.saveSeries(
      {
        title: p.title,
        samples: p.samples,
        ...(p.description !== undefined ? { description: p.description } : {}),
        ...(p.orderingVariable !== undefined ? { ordering_variable: p.orderingVariable } : {}),
      },
      undefined,
      opts,
    );
  },
  onSuccess: (_p, response, qc) => {
    // Splice the new series into the listing and prime the detail cache so the
    // /series/:id navigation has it without a refetch. `isFullSeries` guards
    // against a partial shape (defensive — this mutator only ever resolves
    // with a real POST response, but keep the contract symmetric with
    // saveSeriesMutator's onSuccess).
    if (api.isFullSeries(response)) {
      const existing = qc.getQueryData<Series[]>(queryKeys.seriesList) ?? [];
      qc.setQueryData(queryKeys.seriesList, [...existing, response]);
      qc.setQueryData(queryKeys.series(response.id), response);
    } else {
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
    }
    // Invalidate the scoping proposal sources (same as scopeSeriesMutator).
    qc.invalidateQueries({ queryKey: queryKeys.corpusSampleTags });
    qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
  },
};
