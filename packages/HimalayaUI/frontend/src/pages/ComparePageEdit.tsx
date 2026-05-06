/**
 * ComparePageEdit — edit-mode shell (Plan §Phase 4, Task 4.1 placeholder).
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare/new        — empty draft (create flow)
 *   /experiments/:eid/compare/:id/edit   — edit existing comparison
 *
 * Task 4.1 wires only the routing + sidebar mount. The title/description
 * inputs, draft hydration, and Save / Cancel / Discard buttons land in
 * Task 4.4 once the Zustand draft slot exists (Task 4.3).
 */
import { useParams } from "react-router-dom";
import { ComparisonSidebar } from "../components/ComparisonSidebar";

export function ComparePageEdit(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;
  const scope: "all" | "experiment" = eid !== undefined ? "experiment" : "all";

  return (
    <div
      data-testid="compare-page-edit"
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
      <section className="card overflow-hidden flex-1 min-h-0 flex flex-col items-center justify-center text-fg-muted text-sm p-8 text-center">
        {id === undefined
          ? "Edit mode (new comparison) — title, members, and Save flow land in Task 4.4."
          : `Edit mode for comparison #${id} — title, members, and Save flow land in Task 4.4.`}
      </section>
    </div>
  );
}
