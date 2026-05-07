/**
 * Compare-page URL builders.
 *
 * One source of truth for the route shape used across navigation sites
 * (sidebar pick, edit/fork buttons, conflict modal, lineage badge, forks
 * popover, warm-add menu, needs-review badge). Without this helper the
 * "do I have an :eid in scope?" branching gets duplicated and drifts —
 * for example, the global `/compare/all` scope only had a list URL but no
 * deep-link review URL, so navigating to a comparison from the global
 * sidebar silently fell back to the empty list.
 *
 * URL shape mirrors the routes registered in `AppShell.tsx`:
 *   - `/experiments/:eid/compare`            — experiment-scoped list
 *   - `/experiments/:eid/compare/:id`        — experiment-scoped review
 *   - `/experiments/:eid/compare/:id/edit`   — experiment-scoped edit
 *   - `/experiments/:eid/compare/new`        — experiment-scoped create
 *   - `/compare/all`                         — global list
 *   - `/compare/all/:id`                     — global review
 *   - `/compare/all/:id/edit`                — global edit
 *   - `/compare/all/new`                     — global create
 *
 * The `scope` param is explicit (not derived from `eid`) so callers must
 * decide intentionally whether the URL should pin to an experiment or to
 * the global "all" listing. Mixing them up is a category error a function
 * signature can prevent.
 */

export type CompareScope = "experiment" | "all";

export interface ComparePathOpts {
  scope: CompareScope;
  /** Required for `scope: "experiment"`; ignored for `scope: "all"`. */
  eid?: number | undefined;
  /** Comparison id. Omit for the list/new path. */
  id?: number | undefined;
  /** When true, append `/edit` (only meaningful when `id` is set). */
  edit?: boolean | undefined;
  /** When true, return the create path (`.../new`). Mutually exclusive with `id`. */
  isNew?: boolean | undefined;
}

/**
 * Build a Compare URL from the scope + ids. Throws (TypeScript level) on
 * a missing `eid` for experiment scope; runtime falls back to the global
 * scope so a stray undefined doesn't 404 the user.
 */
export function comparePath(opts: ComparePathOpts): string {
  const { scope, eid, id, edit, isNew } = opts;

  // Build the prefix. When scope=experiment but eid is undefined (defensive
  // — TypeScript can't always prove eid is set), fall through to the global
  // path so the user lands somewhere coherent rather than `/experiments//compare`.
  const prefix =
    scope === "experiment" && eid !== undefined
      ? `/experiments/${eid}/compare`
      : "/compare/all";

  if (isNew) return `${prefix}/new`;
  if (id === undefined) return prefix;
  return edit ? `${prefix}/${id}/edit` : `${prefix}/${id}`;
}
