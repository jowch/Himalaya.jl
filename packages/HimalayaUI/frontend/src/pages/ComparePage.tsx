/**
 * ComparePage — review-mode shell (Plan §Phase 4, Task 4.1).
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare        — sidebar list scoped to experiment, no comparison selected
 *   /experiments/:eid/compare/:id    — review mode of comparison `:id`
 *   /compare/all                     — global listing scope (no experiment context)
 *
 * Layout uses the shared WorkspaceGrid (issue #60) so Compare inherits the
 * same 1400px breakpoint + 700px max-h clamp as IndexPage / InspectPage:
 *   left   = ComparisonSidebar
 *   center = review plot (header + plot + gutter), or empty state
 *   right  = ChatCard (comparison history), or hint when no comparison selected
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useParams, useLocation, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { Skeleton } from "boneyard-js/react";
import { ComparisonSidebar } from "../components/ComparisonSidebar";
import { MultiTracePlot } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { GroupingModeToggle } from "../components/GroupingModeToggle";
import { AnnotationToggles } from "../components/AnnotationToggles";
import { NeedsReviewBadge } from "../components/NeedsReviewBadge";
import { LineageBadge } from "../components/LineageBadge";
import { ForksPopover } from "../components/ForksPopover";
import { ChatCard } from "../components/ChatCard";
import { WorkspaceGrid } from "../components/WorkspaceGrid";
import { FigureExportControls } from "../components/FigureExportControls";
import { HintText } from "../components/ui";
import { buildMultiTraceExportSpec } from "../lib/figure-export/adapters/multiTraceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
import { effectiveGroupingMode } from "../lib/comparison/effectiveGroupingMode";
import {
  useComparison,
  useMemberTraces,
  useMemberTracesLoading,
  useMemberExposures,
  useMemberSamples,
  useExperiment,
} from "../queries";
import { useAppState } from "../state";
import { useCurrentUserId } from "../hooks/useCurrentUserId";
import { comparePath, type CompareScope } from "../lib/comparison/routes";
import { resolveDisplayLabels } from "../lib/comparison/labels";
import type { Comparison, ComparisonMember } from "../api";

// Boneyard fixture for the compare-review-plot skeleton — a stand-in plot
// pane + gutter so the captured bones reflect the dual-column geometry the
// user sees during a true cold fetch (comparison + member-traces both loading).
const COMPARE_PLOT_FIXTURE = (
  <div className="flex-1 min-h-0 flex flex-row gap-3">
    <div className="flex-1 min-w-0 border border-border/40 rounded h-full" />
    <div className="w-[280px] shrink-0 border border-border/40 rounded h-full" />
  </div>
);

export function ComparePage(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const location = useLocation();
  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;
  const scope: "all" | "experiment" =
    location.pathname.startsWith("/compare/all") ? "all" : "experiment";

  return (
    <div
      data-testid="compare-page-body"
      data-scope={scope}
      {...(id !== undefined ? { "data-comparison-id": String(id) } : {})}
      className="flex-1 min-h-0 flex flex-col gap-3 px-4 pb-4 pt-2"
    >
      <WorkspaceGrid
        left={
          <ComparisonSidebar
            experimentId={eid}
            scope={scope}
            activeComparisonId={id}
          />
        }
        center={
          id === undefined ? (
            <div
              data-testid="compare-empty-state"
              className="h-full flex items-center justify-center text-fg-muted text-sm p-8 text-center"
            >
              Pick a comparison from the sidebar, or use{" "}
              <span className="font-medium text-fg ml-1">+ New</span> to create one.
            </div>
          ) : (
            <ReviewPlot id={id} eid={eid} scope={scope} />
          )
        }
        right={
          id === undefined ? (
            <div className="h-full flex items-center justify-center p-4 text-center">
              <HintText>Chat appears once a comparison is selected.</HintText>
            </div>
          ) : (
            <div
              data-testid="compare-review-chat"
              className="h-full min-h-0 flex flex-col"
            >
              <ChatCard entityType="comparison" entityId={id} />
            </div>
          )
        }
        slotClassName={{
          // ComparisonSidebar uses `flex-1`, so the slot needs `display:flex`.
          left:   "flex flex-col min-h-[400px]",
          center: "flex flex-col min-h-[480px]",
          right:  "flex flex-col min-h-[400px]",
        }}
      />
    </div>
  );
}

/**
 * Hosts the multi-trace plot in review mode. Members come from the saved
 * comparison; live `(q, I)` traces are fetched in parallel via
 * `useMemberTraces`.
 */
function ReviewPlot({
  id, eid, scope,
}: {
  id: number;
  eid: number | undefined;
  scope: CompareScope;
}): JSX.Element {
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

  // Phase 9.6 — comparison-level stale flag is the disjunction of per-member
  // is_stale. Hidden when the comparison hasn't loaded yet.
  const isStale = useMemo(() => members.some((m) => m.is_stale), [members]);
  const authorUserId = compQ.data?.created_by ?? null;

  // Phase 9.3 — annotation toggles read straight from Zustand. Both default
  // to `true`; `MultiTracePlot` rebuilds marks when the values change so
  // toggling re-renders without a full plot lifecycle event.
  const showPeakTicks  = useAppState((s) => s.showPeakTicks);
  const showPeakLabels = useAppState((s) => s.showPeakLabels);

  // Phase 9 gap-fix — line-stroke coloring grouping mode + sample-id resolver.
  // The resolver reads from the per-member exposure subscription set up
  // below. Without that explicit subscription the cache never warms in
  // review mode (only the trace key is fetched), and every trace falls back
  // to ORPHAN_FALLBACK gray. See issue #61 + #52.
  //
  // C-4: groupingMode is resolved via effectiveGroupingMode (draft →
  // server record → default). Write goes through setDraftViewGroupingMode
  // which auto-creates an empty draft when none is active (spec §6.4
  // viewer escape hatch).
  const activeDraft = useAppState((s) => s.activeDraft);
  const setDraftViewGroupingMode = useAppState((s) => s.setDraftViewGroupingMode);
  const groupingMode = effectiveGroupingMode(activeDraft, compQ.data);

  // Phase 9.5 — hover/click-to-pin highlight state. `MemberTraceLayer`
  // reads this and recolors that member's confirmed_index peaks to the
  // phase color; non-index peaks stay black. The lifecycle hook below
  // clears the pin on page navigation away (component unmount) and on
  // member-set changes that no longer include the pinned member.
  const highlightedCompareMemberId = useAppState(
    (s) => s.highlightedCompareMemberId,
  );
  const setHighlightedCompareMemberId = useAppState(
    (s) => s.setHighlightedCompareMemberId,
  );
  // Clear stale pin if the pinned member is no longer in the comparison
  // (e.g., re-fetched comparison drops a member). Without this, the pin
  // would persist across navigation between comparisons in the SAME
  // review-mode page.
  useEffect(() => {
    if (highlightedCompareMemberId === undefined) return;
    const stillPresent = members.some((m) => m.id === highlightedCompareMemberId);
    if (!stillPresent) setHighlightedCompareMemberId(undefined);
  }, [
    highlightedCompareMemberId,
    members,
    setHighlightedCompareMemberId,
  ]);
  // Unmount cleanup: navigating away from the compare page (or to edit
  // mode, which is a separate page component) clears the pin.
  useEffect(() => {
    return () => setHighlightedCompareMemberId(undefined);
  }, [setHighlightedCompareMemberId]);

  const exposureIds = useMemo(
    () => members.flatMap((m) => (m.exposure_id !== null ? [m.exposure_id] : [])),
    [members],
  );
  const traces = useMemberTraces(exposureIds);
  const tracesLoading = useMemberTracesLoading(exposureIds);

  // Hydrate per-member exposure rows so `sampleIdFor` and the gutter label
  // resolver below have the data they need (issue #61 + #52). The Map is
  // ref-stable across renders.
  const exposures = useMemberExposures(exposureIds);
  const sampleIds = useMemo(() => {
    const ids = new Set<number>();
    for (const e of exposures.values()) ids.add(e.sample_id);
    return Array.from(ids).sort((a, b) => a - b);
  }, [exposures]);
  const samples = useMemberSamples(sampleIds);

  const sampleIdFor = useCallback(
    (m: ComparisonMember): number | null => {
      if (m.exposure_id === null) return null;
      return exposures.get(m.exposure_id)?.sample_id ?? null;
    },
    [exposures],
  );

  // Resolved per-member labels via the shared resolver (lib/comparison/labels.ts).
  // ComparePageEdit shares the same helper so the fallback chain stays in
  // lockstep across review and edit modes. See issue #69.
  const displayLabelByMemberId = useMemo(
    () => resolveDisplayLabels(members, exposures, samples),
    [members, exposures, samples],
  );

  // Figure export — read experiment name only when in experiment scope.
  const experimentQ = useExperiment(eid ?? 0);
  const experimentName = eid !== undefined
    ? (experimentQ.data?.name ?? `Experiment ${eid}`)
    : undefined;

  const exportFilenameStem = `himalaya-comparison-${
    slugifyForFilename(experimentName ?? "all")
  }-${
    slugifyForFilename(compQ.data?.title ?? "")
  }`;

  const exportSpec = useCallback(() => {
    if (!compQ.data || members.length === 0) {
      throw new Error("FigureExportControls: parent disabled-gate violated");
    }
    return buildMultiTraceExportSpec({
      members,
      traces,
      comparisonTitle: compQ.data.title,
      ...(experimentName !== undefined ? { experimentName } : {}),
      xDomain,
      showPeakTicks,
      showPeakLabels,
      groupingMode,
      sampleIdFor,
      displayLabelByMemberId,
    });
  }, [
    compQ.data, members, traces, experimentName,
    xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
  ]);

  const exportDisabled =
    members.length === 0
    || traces.size === 0
    || members.every((m) => m.exposure_id === null);

  // Cold-fetch gate. Per CLAUDE.md boneyard rules: gate on `query.isLoading`
  // not `isPending` so disabled queries / background refetches don't flicker.
  const plotLoading = compQ.isLoading || tracesLoading;

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
  }, [plotLoading]);

  return (
    <div className="flex-1 min-h-0 flex flex-col p-4 gap-3" data-testid="compare-review-plot">
      <div
        data-testid="compare-review-header"
        className="flex items-center gap-3 flex-wrap"
      >
        <GroupingModeToggle mode={groupingMode} onChange={setDraftViewGroupingMode} />
        <AnnotationToggles />
        {isStale && (
          <NeedsReviewBadge
            comparisonId={id}
            experimentId={eid}
            scope={scope}
            authorUserId={authorUserId}
          />
        )}
        {compQ.data && (
          <EditOrForkButton
            comparison={compQ.data}
            experimentId={eid}
            scope={scope}
          />
        )}
        {compQ.data && (
          <LineageBadge comparison={compQ.data} experimentId={eid} scope={scope} />
        )}
        <ForksPopover comparisonId={id} experimentId={eid} scope={scope} />
        <span className="ml-auto inline-flex items-center gap-1">
          <FigureExportControls
            spec={exportSpec}
            filenameStem={exportFilenameStem}
            ariaContext="comparison plot"
            disabled={exportDisabled}
          />
        </span>
      </div>
      <Skeleton
        name="compare-review-plot"
        className="flex-1 min-h-0 flex flex-row gap-3"
        loading={plotLoading}
        stagger={50}
        transition={200}
        fixture={COMPARE_PLOT_FIXTURE}
        fallback={<div className="flex-1 flex items-center justify-center"><HintText>Loading comparison…</HintText></div>}
      >
        <div ref={plotColRef} className="flex-1 min-w-0">
          <MultiTracePlot
            members={members}
            traces={traces}
            xDomain={xDomain}
            onXDomain={setXDomain}
            showPeakTicks={showPeakTicks}
            showPeakLabels={showPeakLabels}
            highlightedMemberId={highlightedCompareMemberId}
            groupingMode={groupingMode}
            sampleIdFor={sampleIdFor}
          />
        </div>
        <div
          className="w-[280px] shrink-0 relative"
          data-testid="compare-review-gutter"
        >
          <MemberMetaGutter
            members={members}
            panelHeight={panelHeight}
            mode="review"
            displayLabelByMemberId={displayLabelByMemberId}
          />
        </div>
      </Skeleton>
    </div>
  );
}

/**
 * Author-vs-fork affordance (Plan §Phase 11, Task 11.2). Mutually exclusive
 * — Edit when the current user authored the comparison, Fork otherwise. The
 * orphan-author case (`comparison.created_by === null`) shows Fork to ALL
 * users since no one matches null; combined with the backend's spec
 * §Authorship "fork-only" gate, this is the right fallback.
 *
 * Edit click navigates to the edit-mode shell and seeds the Zustand draft
 * via `loadDraftFromComparison` so the editor has the full saved state to
 * mutate. Fork click creates a brand-new draft pre-populated from the
 * parent's data + lineage (`forkedFromId` + `forkedAtHash`) and navigates
 * to the create flow; submit will POST /api/comparisons with the lineage
 * fields per Phase 3's `SaveComparisonBody` contract.
 */
function EditOrForkButton({
  comparison, experimentId, scope,
}: {
  comparison: Comparison;
  experimentId: number | undefined;
  scope: CompareScope;
}): JSX.Element {
  const navigate = useNavigate();
  const qc = useQueryClient();
  const currentUserId = useCurrentUserId();
  const loadDraft = useAppState((s) => s.loadDraftFromComparison);
  const startFork = useAppState((s) => s.startForkDraft);

  const isAuthor =
    comparison.created_by !== null
    && currentUserId !== undefined
    && currentUserId === comparison.created_by;

  if (isAuthor) {
    const onEdit = (): void => {
      loadDraft(comparison, qc);
      navigate(comparePath({ scope, eid: experimentId, id: comparison.id }));
    };
    return (
      <button
        type="button"
        data-testid="comparison-edit"
        onClick={onEdit}
        className="px-1.5 py-0.5 rounded text-xs text-fg-dim hover:text-fg
                   hover:bg-bg-hover border border-transparent hover:border-border"
      >
        Edit
      </button>
    );
  }

  const onFork = (): void => {
    startFork(comparison, qc);
    navigate(comparePath({ scope, eid: experimentId, isNew: true }));
  };
  return (
    <button
      type="button"
      data-testid="comparison-fork"
      onClick={onFork}
      className="px-1.5 py-0.5 rounded text-xs text-fg-dim hover:text-fg
                 hover:bg-bg-hover border border-transparent hover:border-border"
    >
      Fork
    </button>
  );
}
