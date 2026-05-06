/**
 * ComparePage — review-mode shell (Plan §Phase 4, Task 4.1).
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare        — sidebar list scoped to experiment, no comparison selected
 *   /experiments/:eid/compare/:id    — review mode of comparison `:id`
 *   /compare/all                     — global listing scope (no experiment context)
 *
 * Layout mirrors the three-card workspace pattern (sidebar | main).
 * The plot, member panel, chat, badges, and edit/fork affordances are
 * built out across Phases 6–11; this file is only the shell that hosts them.
 */
import { useParams, useLocation } from "react-router-dom";
import { ComparisonSidebar } from "../components/ComparisonSidebar";

export function ComparePage(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const location = useLocation();
  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;
  const scope: "all" | "experiment" =
    location.pathname.startsWith("/compare/all") ? "all" : "experiment";

  return (
    <div
      data-testid="compare-page"
      data-scope={scope}
      {...(id !== undefined ? { "data-comparison-id": String(id) } : {})}
      className="flex-1 min-h-0 flex gap-3 px-4 pb-4 pt-2"
    >
      <aside className="card overflow-hidden w-[300px] shrink-0 flex flex-col">
        <ComparisonSidebar
          experimentId={eid}
          scope={scope}
          activeComparisonId={id}
        />
      </aside>
      <section className="card overflow-hidden flex-1 min-h-0 flex flex-col">
        {id === undefined ? (
          <div
            data-testid="compare-empty-state"
            className="flex-1 flex items-center justify-center text-fg-muted text-sm p-8 text-center"
          >
            Pick a comparison from the sidebar, or use{" "}
            <span className="font-medium text-fg ml-1">+ New</span> to create one.
          </div>
        ) : (
          <div
            data-testid="compare-review-placeholder"
            className="flex-1 flex items-center justify-center text-fg-muted text-sm p-8 text-center"
          >
            Comparison #{id} (review mode) — plot + chat + badges land in later phases.
          </div>
        )}
      </section>
    </div>
  );
}
