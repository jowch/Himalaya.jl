/**
 * commitSeriesPlate mutator (#167). Commits the plate via
 * `POST /api/series/:id/commit` — the old "submit". No optimistic write
 * (spinner). `series_plate_committed` is the one series event carrying a
 * `post_state` envelope (the full `fetch_series_with_plate` projection), so
 * `synthesizeFromSse` can return a complete Series on the SSE-wins path.
 *
 * A 409 (content_hash conflict) surfaces as a generic `ApiError`; the typed
 * conflict modal is I3.5b's concern.
 */
import * as api from "../../../api";
import type {
  AuthOpts, Series, SeriesMemberInput, CommitSeriesPlateBody,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface CommitSeriesPlateInput {
  id: number;
  members: SeriesMemberInput[];
  expected_content_hash?: string;
}

interface CommitSeriesPlateScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const commitSeriesPlateMutator: Mutator<CommitSeriesPlateInput, CommitSeriesPlateScope, Series> = {
  kind: "series_commit",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => {
    // Annotate with the shared CommitSeriesPlateBody type — do not re-declare
    // it inline (drift risk).
    const body: CommitSeriesPlateBody = { members: p.members };
    if (p.expected_content_hash !== undefined) {
      body.expected_content_hash = p.expected_content_hash;
    }
    return api.commitSeriesPlate(p.id, body, buildAuthOpts(p));
  },
  onSuccess: (_p, response, qc) => {
    // `api.isFullSeries` distinguishes a real HTTP Series from the partial
    // SSE-wins shape; fall back to invalidate when the shape is incomplete.
    // The else branch reads `id` via an explicit cast — the guard narrows
    // `response` to `never` there.
    if (api.isFullSeries(response)) {
      qc.setQueryData(queryKeys.series(response.id), response);
    } else {
      const id = (response as { id?: unknown }).id;
      if (typeof id === "number") {
        qc.invalidateQueries({ queryKey: queryKeys.series(id) });
      }
    }
    qc.invalidateQueries({ queryKey: queryKeys.seriesList });
  },
  synthesizeFromSse: (remote, _base) => {
    // post_state IS the full fetch_series_with_plate projection (master plan
    // §5.2). Return it RAW — do NOT spread `_base` (QueueResponseMeta) in:
    // `onSuccess`'s looksFull branch writes this object straight to the detail
    // cache via setQueryData, and merging event_id / client_op_id /
    // analysis_inputs_hash would pollute the cached Series row and diverge
    // from the clean post_state splice in the applyRemoteToCache arm. `_base`
    // is intentionally unused (underscore-prefixed for noUnusedParameters).
    if (remote.post_state != null) {
      return remote.post_state as Series;
    }
    return undefined;
  },
};
