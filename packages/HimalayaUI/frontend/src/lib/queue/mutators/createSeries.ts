import type { Mutator, RollbackContext } from "../types";
import type { AuthOpts, Series, SeriesSampleInput } from "../../../api";
import { authOpts } from "../../authOpts";
import { queryKeys } from "../../../queries";
import * as api from "../../../api";

export interface CreateSeriesInput {
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  ordering_variable?: string | null;
}

interface CreateSeriesScope {
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
 * createSeriesMutator — series CREATE only (M-A: Op B of the scoping
 * scope→create chain). A SINGLE write: `api.saveSeries(body, undefined, …)`
 * (no id → POST /api/series).
 *
 * WHY A SINGLE-WRITE MUTATOR. The queue resolves an op's deferred (keyed on the
 * op's clientOpId) with WHICHEVER lands first — the HTTP `request` return OR an
 * own-op SSE frame on that same op-id (see useQueueMutation's mutationFn). A
 * compound two-write op (the abandoned scopeAndCreateSeries) emitted `add_tag`
 * frames from its first write under the bare op-id; those resolved the deferred
 * BEFORE the in-flight create returned, so `mutation.data` was the tag
 * confirmation, not the Series. Splitting create into its own single-write op
 * means its ONLY own-op SSE frame is `series_created` (matching this op's
 * op-id), so the deferred reliably resolves with the full Series —
 * `mutation.data` is the created Series (exactly like saveSeriesMutator, which
 * the builder reads `.data.members` from). The tag write is a SEPARATE queue
 * op (scopeSeriesMutator), sequenced ahead of this one by the page's ref state
 * machine.
 *
 * kind='series_save' is the frontend OpKind (the create route emits a
 * `series_created` wire frame). onMutate is a no-op (the user navigates away on
 * success). No `synthesizeFromSse` — foreign `series_created` frames route
 * through `saveSeriesMutator.synthesizeFromSse` via mutatorRegistry (the
 * registry never resolves a foreign event to this mutator).
 *
 * onSuccess: splices the new Series into the listing cache and primes the
 * detail cache (so the navigation to /series/:id reads it immediately), then
 * invalidates the two scoping proposal sources (the corpus tag + picker
 * projections, which the just-completed tag write changed).
 */
export const createSeriesMutator: Mutator<CreateSeriesInput, CreateSeriesScope, Series> = {
  kind: "series_save",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) =>
    api.saveSeries(
      {
        title: p.title,
        samples: p.samples,
        ...(p.description !== undefined ? { description: p.description } : {}),
        ...(p.ordering_variable !== undefined ? { ordering_variable: p.ordering_variable } : {}),
      },
      undefined, // no id → POST /api/series (create)
      buildAuthOpts(p),
    ),
  onSuccess: (_p, response, qc) => {
    // Splice the new series into the listing and prime the detail cache so the
    // /series/:id navigation has it without a refetch. `isFullSeries` guards
    // against a partial shape (defensive; this single-write op resolves with a
    // real POST response, but keep the contract symmetric with
    // saveSeriesMutator's onSuccess).
    if (api.isFullSeries(response)) {
      const existing = qc.getQueryData<Series[]>(queryKeys.seriesList) ?? [];
      qc.setQueryData(queryKeys.seriesList, [...existing, response]);
      qc.setQueryData(queryKeys.series(response.id), response);
    } else {
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
    }
    // Invalidate the scoping proposal sources (same as scopeSeriesMutator) —
    // the tag write that preceded this create changed them.
    qc.invalidateQueries({ queryKey: queryKeys.corpusSampleTags });
    qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
  },
};
