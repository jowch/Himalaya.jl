/**
 * Series builder draft — type definitions and sessionStorage persistence
 * (I3.5b). A separate namespace from the comparison `ActiveDraft`: a series
 * recipe edits `series_samples` MEMBERSHIP (sample_id, position, pinned,
 * excluded), a different shape from a comparison's per-member plate. Holding
 * them in one slot would let the two flows clobber each other across tabs.
 *
 * The commit route is last-write-wins (Plan 6a): the draft no longer carries a
 * `baseHash` / `expected_content_hash` — the backend stopped 409ing on a stale
 * hash, so nothing reads it.
 *
 * Optionality is encoded as `T | null` / `T | undefined` (not `T?`) per the
 * `exactOptionalPropertyTypes` gotcha so round-trip Zustand `set` works.
 *
 * **Module structure note:** like `comparison/draft.ts`, this file is loaded by
 * `state.ts` at module-init (via `loadSeriesDraftFromSession()`), so it MUST
 * stay free of any transitive `queries.ts` import.
 */
import type { OrderRule } from "../../api";

/**
 * One editable membership row. `id` is the `series_samples.id` when persisted
 * (positive) or a `nextOptimisticId()` negative placeholder when freshly added
 * and not yet round-tripped. The id is a LOCAL render-key only — never sent to
 * the server (the wire `SeriesSampleInput` is positional, id-less) and never
 * used to reconcile state (replay-volatile, master plan §11; reconcile by
 * `(series_id, position)`).
 */
export interface SeriesRecipeRow {
  id: number;
  sample_id: number;
  position: number;
  pinned: boolean;
  excluded: boolean;
}

export interface SeriesDraft {
  /** The existing series being edited. Always a real positive id — there is
   *  no "create a series from the builder" path (creation is scoping, I3.4). */
  id: number;
  title: string;
  description: string;
  orderingVariable: string | null;
  orderRule: OrderRule;
  recipe: SeriesRecipeRow[];
}

export type SeriesDraftSlot = SeriesDraft | null;

export const SERIES_DRAFT_KEY = "himalaya-ui:series-draft";
export const SERIES_DRAFT_VERSION = 1;

interface PersistedEnvelope {
  version: number;
  draft: SeriesDraft;
}

export function emptySeriesDraft(id: number): SeriesDraft {
  return {
    id,
    title: "",
    description: "",
    orderingVariable: null,
    orderRule: "manual",
    recipe: [],
  };
}

export function persistSeriesDraftToSession(draft: SeriesDraft | null): void {
  try {
    if (draft === null) {
      sessionStorage.removeItem(SERIES_DRAFT_KEY);
      return;
    }
    const env: PersistedEnvelope = { version: SERIES_DRAFT_VERSION, draft };
    sessionStorage.setItem(SERIES_DRAFT_KEY, JSON.stringify(env));
  } catch {
    // sessionStorage can throw under quota/access restrictions; best-effort.
  }
}

export function loadSeriesDraftFromSession(): SeriesDraft | null {
  try {
    const raw = sessionStorage.getItem(SERIES_DRAFT_KEY);
    if (!raw) return null;
    const parsed = JSON.parse(raw) as Partial<PersistedEnvelope> | null;
    if (!parsed || typeof parsed !== "object") return null;
    if (parsed.version !== SERIES_DRAFT_VERSION) return null;
    if (!parsed.draft || typeof parsed.draft !== "object") return null;
    return parsed.draft as SeriesDraft;
  } catch {
    return null;
  }
}
