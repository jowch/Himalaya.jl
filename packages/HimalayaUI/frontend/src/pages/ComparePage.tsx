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
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useParams, useLocation } from "react-router-dom";
import { ComparisonSidebar } from "../components/ComparisonSidebar";
import { MultiTracePlot } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { GroupingModeToggle } from "../components/GroupingModeToggle";
import { useComparison, useMemberTraces } from "../queries";
import { useAppState } from "../state";

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
          <ReviewPlot id={id} />
        )}
      </section>
    </div>
  );
}

/**
 * Hosts the multi-trace plot in review mode. Members come from the saved
 * comparison; live `(q, I)` traces are fetched in parallel via
 * `useMemberTraces`.
 */
function ReviewPlot({ id }: { id: number }): JSX.Element {
  const compQ = useComparison(id);
  // Per-comparison zoom keying — selecting only the slice for `id` so this
  // component does not re-render on zoom changes to other comparisons.
  const xDomain = useAppState((s) => s.compareXDomains[id] ?? null);
  const setCompareXDomain = useAppState((s) => s.setCompareXDomain);
  const setXDomain = useCallback(
    (d: [number, number] | null) => setCompareXDomain(id, d),
    [setCompareXDomain, id],
  );

  const members = useMemo(() => {
    if (!compQ.data) return [];
    return [...compQ.data.members].sort((a, b) => a.display_order - b.display_order);
  }, [compQ.data]);

  const exposureIds = useMemo(
    () => members.flatMap((m) => (m.exposure_id !== null ? [m.exposure_id] : [])),
    [members],
  );
  const traces = useMemberTraces(exposureIds);

  // Track the plot column's height so the gutter rows align pixel-for-pixel
  // with the plot's y-bands (both consumers share `computeYBands`).
  const plotColRef = useRef<HTMLDivElement>(null);
  const [panelHeight, setPanelHeight] = useState(0);
  useEffect(() => {
    const el = plotColRef.current;
    if (!el) return;
    setPanelHeight(el.clientHeight);
    const obs = new ResizeObserver(() => {
      if (plotColRef.current) setPanelHeight(plotColRef.current.clientHeight);
    });
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

  return (
    <div className="flex-1 min-h-0 flex flex-col p-4 gap-3" data-testid="compare-review-plot">
      <div
        data-testid="compare-review-header"
        className="flex items-center gap-2 flex-wrap"
      >
        <GroupingModeToggle />
      </div>
      <div className="flex-1 min-h-0 flex flex-row gap-3">
        <div ref={plotColRef} className="flex-1 min-w-0">
          <MultiTracePlot
            members={members}
            traces={traces}
            xDomain={xDomain}
            onXDomain={setXDomain}
          />
        </div>
        <div
          className="w-[280px] shrink-0 relative"
          data-testid="compare-review-gutter"
        >
          <MemberMetaGutter members={members} panelHeight={panelHeight} mode="review" />
        </div>
      </div>
    </div>
  );
}
