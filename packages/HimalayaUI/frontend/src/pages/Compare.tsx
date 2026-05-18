/**
 * Compare.tsx — unified single-mode shell for all compare routes.
 *
 * Compare UX C-15: this file folds the former `ComparePage` (review) and
 * `ComparePageEdit` (create / edit) bodies into one. The top-level `Compare`
 * component picks which body to mount:
 *
 *   - URL ends `/new`  ⇒ edit body (the hydration effect mints a fresh draft).
 *   - a draft is active ⇒ edit body (inline edit gesture on `/compare/:id`).
 *   - otherwise         ⇒ review body.
 *
 * The two bodies keep separate hook graphs on purpose — review reads the
 * server `Comparison`, edit reads the Zustand `activeDraft` — so they live as
 * sibling components inside this file rather than one branching subtree. The
 * shared header components (`CompareTitleStrip` / `CompareStatusSurface` /
 * `CompareToolbar` / `SavePill`) and the shared plot/gutter mounts are used by
 * both. Structural testids are preserved verbatim from the old two-file split
 * so the pre-C-15 test suite stays green.
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare        — review, no comparison selected
 *   /experiments/:eid/compare/:id    — review of comparison `:id` (or edit
 *                                       once a draft is active)
 *   /experiments/:eid/compare/new    — create flow (empty draft)
 *   /compare/all[...]                — global listing scope
 *
 * The wrapper div carries `data-testid="compare-page"` (the canonical testid
 * for the compare surface) and `className="contents"` so it is transparent to
 * the flex/grid layout of its parent — it does not consume a grid track in
 * WorkspaceGrid's 3-column subgrid.
 */
import {
  useCallback, useEffect, useMemo, useRef, useState,
} from "react";
import { useParams, useLocation, useNavigate } from "react-router-dom";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { Skeleton } from "boneyard-js/react";
import { ComparisonSidebar } from "../components/ComparisonSidebar";
import { ComparisonPickerPanel } from "../components/ComparisonPickerPanel";
import { MultiTracePlot } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { ActiveBandProvider } from "../components/ActiveBandContext";
import { GroupingModeToggle } from "../components/GroupingModeToggle";
import { AnnotationToggles } from "../components/AnnotationToggles";
import { CompareTitleStrip } from "../components/CompareTitleStrip";
import { CompareStatusSurface } from "../components/CompareStatusSurface";
import { CompareToolbar } from "../components/CompareToolbar";
import { SavePill } from "../components/SavePill";
import { ChatCard } from "../components/ChatCard";
import { WorkspaceGrid } from "../components/WorkspaceGrid";
import { FigureExportControls } from "../components/FigureExportControls";
import { HintText } from "../components/ui";
import { buildMultiTraceExportSpec } from "../lib/figure-export/adapters/multiTraceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
import { effectiveGroupingMode } from "../lib/comparison/effectiveGroupingMode";
import { comparePath, type CompareScope } from "../lib/comparison/routes";
import { resolveDisplayLabels } from "../lib/comparison/labels";
import { computeMemberSnapshot } from "../lib/comparison/snapshot";
import { prefetchColdMembers } from "../lib/comparison/prefetchMembers";
import {
  useComparison,
  useComparisonForks,
  useMemberTraces,
  useMemberTracesLoading,
  useMemberExposures,
  useMemberSamples,
  useExperiment,
  useDeleteComparison,
  useSaveComparison,
} from "../queries";
import { useAppState } from "../state";
import { useCurrentUserId } from "../hooks/useCurrentUserId";
import { useCompareMode } from "../hooks/useCompareMode";
import * as api from "../api";
import type {
  Comparison, ComparisonMemberInput, SaveComparisonBody, SeriesMember,
} from "../api";
import type { DraftMember } from "../lib/comparison/draft";

// ── Top-level shell ─────────────────────────────────────────────────────────

export function Compare(): JSX.Element {
  const location = useLocation();
  const isNew = location.pathname.endsWith("/new");
  const hasDraft = useAppState((s) => s.activeDraft !== null);
  // `/new` always lands on the edit body (its hydration effect mints a draft);
  // an active draft on `/compare/:id` also routes to the edit body (the inline
  // edit gesture). Everything else is review.
  const edit = isNew || hasDraft;
  return (
    <div data-testid="compare-page" className="contents">
      {edit ? <EditBody /> : <ReviewBody />}
    </div>
  );
}

// ── useStaleAgainstHash (Compare UX C-15 Step 4b) ───────────────────────────
//
// SSE → viewing-stale detection. Tracks the comparison `content_hash` the
// user last "saw"; while a foreign update drifts it, `previousHash` is the
// hash the banner should warn against. Three behavioral pins:
//   (1) first sighting anchors `acked` — no banner flash on cold load;
//   (2) `acknowledge()` rebases `acked` → banner clears;
//   (3) own-op completion (draft non-null → null) rebases `acked` so the
//       save's own content_hash bump does not flash the banner.
function useStaleAgainstHash(comparison: Comparison | undefined): {
  previousHash: string | undefined;
  acknowledge: () => void;
} {
  const draft = useAppState((s) => s.activeDraft);
  const [acked, setAcked] = useState<string | undefined>(undefined);

  useEffect(() => {
    if (comparison === undefined) return;
    if (acked === undefined) {
      setAcked(comparison.content_hash);
      return;
    }
    // No auto-advancement: acked advances only via acknowledge() or
    // own-op-rebase.
  }, [comparison?.content_hash, acked]);

  // Pin (3) — own-op no-flash — is delivered by the remount path, not an
  // effect. `ReviewBody` (and with it this hook) unmounts the instant a
  // draft activates and remounts once the draft clears, so the
  // first-sighting branch above re-anchors `acked` to the post-save hash.
  // An explicit non-null → null draft-rebase effect was dropped (Compare UX
  // review): under the ReviewBody/EditBody split it could never fire
  // in-place. Reinstate it if this hook is ever hoisted above that split.

  const acknowledge = useCallback(() => {
    if (comparison !== undefined) setAcked(comparison.content_hash);
  }, [comparison?.content_hash]);

  if (draft !== null) return { previousHash: undefined, acknowledge };
  if (comparison === undefined) return { previousHash: undefined, acknowledge };
  if (acked === undefined) return { previousHash: undefined, acknowledge };
  return {
    previousHash: acked === comparison.content_hash ? undefined : acked,
    acknowledge,
  };
}

// ── Review body (was ComparePage) ───────────────────────────────────────────

// Boneyard fixture for the compare-review-plot skeleton — a stand-in plot
// pane + gutter so the captured bones reflect the dual-column geometry the
// user sees during a true cold fetch (comparison + member-traces both loading).
const COMPARE_PLOT_FIXTURE = (
  <div className="flex-1 min-h-0 flex flex-row gap-3">
    <div className="flex-1 min-w-0 border border-border/40 rounded h-full" />
    <div className="w-[280px] shrink-0 border border-border/40 rounded h-full" />
  </div>
);

function ReviewBody(): JSX.Element {
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

  const members = useMemo<SeriesMember[]>(() => {
    if (!compQ.data) return [];
    // I3.2 bridge — comparison members lack series_id; the render pipeline
    // never reads it. Deleted with this page at I3.6.
    return [...compQ.data.members]
      .sort((a, b) => a.display_order - b.display_order)
      .map((m) => ({ ...m, series_id: 0 }));
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
  // mode, which swaps in the edit body) clears the pin.
  useEffect(() => {
    return () => setHighlightedCompareMemberId(undefined);
  }, [setHighlightedCompareMemberId]);

  // ── Compare UX C-12 — review-mode header wiring ──────────────────────────
  // Author byline + isCurrentUserAuthor. `Comparison` carries only the
  // numeric `created_by`; resolve the display name from the shared `["users"]`
  // cache (same source `useCurrentUserId` reads). No dedicated `useUser`
  // hook exists, so the lookup is inlined here.
  const navigate = useNavigate();
  const qc = useQueryClient();
  const currentUserId = useCurrentUserId();
  const usersQ = useQuery({
    queryKey: ["users"] as const,
    queryFn: () => api.listUsers(),
  });
  const authorUsername = useMemo(() => {
    if (authorUserId === null) return null;
    return (usersQ.data ?? []).find((u) => u.id === authorUserId)?.username
      ?? null;
  }, [usersQ.data, authorUserId]);
  const isCurrentUserAuthor =
    authorUserId !== null
    && currentUserId !== undefined
    && currentUserId === authorUserId;

  // Fork lineage — surface "forked from <title>" with a deep link to the
  // parent's review page when the parent still exists.
  const forkedFromId = compQ.data?.forked_from_id ?? null;
  const forkedFromHref = forkedFromId !== null
    ? comparePath({ scope, eid, id: forkedFromId })
    : null;

  // Toolbar action handlers. `EditOrForkButton` / `NeedsReviewBadge` /
  // `ForksPopover` callsites are gone (Compare UX C-12); the toolbar's
  // overflow menu + status surface now own these gestures.
  const forksQ = useComparisonForks(id);
  // Pre-resolve each child fork to its review-page path so CompareToolbar
  // stays a dumb leaf (Compare UX C-17 — the old ForksPopover's job).
  const forksList = useMemo(
    () => (forksQ.data ?? []).map((f) => ({
      id: f.id,
      title: f.title,
      href: comparePath({ scope, eid, id: f.id }),
    })),
    [forksQ.data, scope, eid],
  );
  const startFork = useAppState((s) => s.startForkDraft);
  const deleteMut = useDeleteComparison();

  const handleCopyLink = useCallback(() => {
    void navigator.clipboard?.writeText(window.location.href);
  }, []);

  const handleDelete = useCallback(() => {
    deleteMut.mutate({ id });
    navigate(comparePath({ scope, eid }));
  }, [deleteMut, id, navigate, scope, eid]);

  // Fork (was EditOrForkButton's non-author branch): seed a brand-new draft
  // carrying the parent's lineage and navigate to the create flow.
  const handleFork = useCallback(() => {
    if (!compQ.data) return;
    startFork(compQ.data, qc);
    navigate(comparePath({ scope, eid, isNew: true }));
  }, [compQ.data, startFork, qc, navigate, scope, eid]);

  // Re-snapshot (was NeedsReviewBadge): the author navigates to the bare
  // comparison URL — Compare UX Phase B made that the inline edit surface,
  // where `loadDraftFromComparison` recomputes snapshots against the cache.
  const handleReanalyze = useCallback(() => {
    navigate(comparePath({ scope, eid, id }));
  }, [navigate, scope, eid, id]);

  // ── Compare UX C-15 Step 4b — SSE → viewing-stale detection ───────────────
  const { previousHash, acknowledge } = useStaleAgainstHash(compQ.data);
  // `useCompareMode` derives the mode from draft + comparison + the stale
  // hash. Review mode never reaches editing-*; the hook is wired here so the
  // viewing-stale tag drives the server-update banner below.
  useCompareMode({
    comparison: compQ.data,
    currentUserId,
    staleAgainstHash: previousHash,
  });

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
    (m: SeriesMember): number | null => {
      if (m.exposure_id === null) return null;
      return exposures.get(m.exposure_id)?.sample_id ?? null;
    },
    [exposures],
  );

  // Resolved per-member labels via the shared resolver (lib/comparison/labels.ts).
  // The edit body shares the same helper so the fallback chain stays in
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
        className="flex flex-col gap-2"
      >
        <CompareTitleStrip
          title={compQ.data?.title ?? ""}
          description={compQ.data?.description ?? null}
          memberCount={members.length}
          authorUsername={authorUsername}
          isCurrentUserAuthor={isCurrentUserAuthor}
          lastEventAt={compQ.data?.last_event_at ?? null}
          forkedFromTitle={compQ.data?.forked_from_title ?? null}
          forkedFromHref={forkedFromHref}
          onTitleChange={() => {}}
          onDescChange={() => {}}
          readOnly
        />
        <CompareStatusSurface
          needsReview={isStale ? { onReanalyze: handleReanalyze } : null}
          serverUpdate={
            previousHash !== undefined
              ? { previousHash, onAcknowledge: acknowledge }
              : null
          }
          savedAt={null}
        />
        <CompareToolbar
          groupingControl={
            <GroupingModeToggle
              mode={groupingMode}
              onChange={setDraftViewGroupingMode}
            />
          }
          annotationControl={<AnnotationToggles />}
          forksList={forksList}
          onCopyLink={handleCopyLink}
          onDelete={handleDelete}
          onDiscardChanges={null}
          onFork={handleFork}
          exportControl={
            <FigureExportControls
              spec={exportSpec}
              filenameStem={exportFilenameStem}
              ariaContext="comparison plot"
              disabled={exportDisabled}
            />
          }
          saveControl={null}
        />
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

// ── Edit body (was ComparePageEdit) ─────────────────────────────────────────

/**
 * Convert a draft member into a SeriesMember-shaped object suitable for
 * `MultiTracePlot`. Unsaved drafts have `id = undefined`; we substitute a
 * stable negative synthetic id keyed by display_order so the plot's per-member
 * keying stays consistent across re-renders. Snapshot can also be undefined
 * mid-edit; the plot tolerates a null snapshot (no peaks rendered).
 */
function draftToMember(d: DraftMember): SeriesMember {
  return {
    id: d.id ?? -(d.display_order + 1),
    series_id: 0,
    exposure_id: d.exposure_id,
    display_order: d.display_order,
    band_height: d.band_height,
    y_offset: d.y_offset,
    normalization: d.normalization,
    color_override: d.color_override ?? null,
    label_override: d.label_override ?? null,
    q_window_min: d.q_window_min ?? null,
    q_window_max: d.q_window_max ?? null,
    peak_display: d.peak_display ?? null,
    snapshot: d.snapshot ?? null,
    is_stale: false,
    created_by: null,
    created_at: null,
  };
}

function EditBody(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const navigate = useNavigate();
  const location = useLocation();
  const qc = useQueryClient();

  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;
  // Pathname-derived scope mirrors the review body. The picker / sidebar /
  // post-save navigation all branch on this so /compare/all/:id stays in the
  // global scope rather than silently jumping back to /compare/all.
  const scope: "all" | "experiment" =
    location.pathname.startsWith("/compare/all") ? "all" : "experiment";

  const draft = useAppState((s) => s.activeDraft);
  const startNewDraft       = useAppState((s) => s.startNewDraft);
  const loadDraftFromComp   = useAppState((s) => s.loadDraftFromComparison);
  const setDraftTitle       = useAppState((s) => s.setDraftTitle);
  const setDraftDescription = useAppState((s) => s.setDraftDescription);
  const discardDraft        = useAppState((s) => s.discardDraft);
  const setHighlightedCompareMemberId = useAppState(
    (s) => s.setHighlightedCompareMemberId,
  );

  // Phase 9.5 — entering edit mode always clears any review-mode pin so
  // peak-click cycle isn't confused by a hovered/pinned highlight from
  // a sibling tab's review session. Mirrors the same lifecycle as
  // `discardDraft` on Cancel/Save.
  useEffect(() => {
    setHighlightedCompareMemberId(undefined);
  }, [setHighlightedCompareMemberId]);

  // Hydrate draft from URL on mount / when URL id changes.
  //
  // Both actions are idempotent (see state.ts):
  //   - `startNewDraft` is a no-op when the active draft is already a
  //     fresh one (id undefined) — preserves an in-progress create.
  //   - `loadDraftFromComp` is a no-op when the active draft already
  //     matches the comparison id — preserves an in-progress edit.
  // The guards live in the action bodies, so this effect can stay
  // dep-exhaustive without reading `draft` and without an eslint-disable.
  const comparisonQ = useComparison(id);
  useEffect(() => {
    if (id === undefined) {
      // /new: start an empty draft (no-op if already on a fresh draft).
      startNewDraft();
    } else if (comparisonQ.data) {
      // /:id: load from the server fetch with cache-derived snapshots
      // (no-op if the draft already matches this comparison id).
      loadDraftFromComp(comparisonQ.data, qc);
    }
  }, [id, comparisonQ.data, startNewDraft, loadDraftFromComp, qc]);

  const pickerSearchRef = useRef<HTMLInputElement>(null);

  const save = useSaveComparison();
  const pendingSubmitRef = useRef(false);

  // ── Compare UX C-13/C-14 — `CompareMode` drives the SavePill copy variant
  // and the fork-morph save branch. Declared above `handleSave` because that
  // handler branches on `compareMode.kind === "editing-as-fork-of"`.
  const currentUserId = useCurrentUserId();
  const compareMode = useCompareMode({
    comparison: comparisonQ.data,
    currentUserId,
  });

  const goToReview = useCallback(
    (newId: number) => {
      navigate(comparePath({ scope, eid, id: newId }));
    },
    [navigate, scope, eid],
  );

  const goToList = useCallback(() => {
    navigate(comparePath({ scope, eid }));
  }, [navigate, scope, eid]);

  const handleCancel = useCallback(() => {
    if (id !== undefined) {
      // Cancel from edit-existing → return to that comparison's review page.
      navigate(comparePath({ scope, eid, id }));
    } else {
      // Cancel from create → return to list.
      goToList();
    }
  }, [navigate, id, scope, eid, goToList]);

  // Compare UX C-13 — "Discard changes" now lives inside the toolbar's ⋯-more
  // menu (a less-prominent gesture than the old standalone button). Guard it
  // with a confirm so a stray menu click can't silently throw away the draft.
  const handleDiscard = useCallback(() => {
    if (!window.confirm("Discard all unsaved changes to this comparison?")) return;
    discardDraft();
    goToList();
  }, [discardDraft, goToList]);

  // Build the save payload fresh from a draft snapshot. Kept as a pure
  // function (not a closure over `draft`) so the C-14 fork flow can call it
  // with the *morphed* draft re-read from the store — see handleSave.
  //
  // CRITICAL (idempotency): each `save.mutate()` call must receive its own
  // freshly-built payload object. `useQueueMutation.mutationFn` mints a
  // `client_op_id` per `mutate()` call; reusing a single stale payload would
  // route the post-morph submission through the pre-morph `id` and the
  // create-path route would never fire.
  const buildSavePayload = useCallback(
    (d: typeof draft): (SaveComparisonBody & { id?: number }) | null => {
      if (d === null) return null;
      // Compute a fresh snapshot per member at submit time (Plan §Task 4.3).
      const members: ComparisonMemberInput[] = d.members.map((m) => {
        const snapshot = m.exposure_id !== null
          ? computeMemberSnapshot(m.exposure_id, qc)
          : (m.snapshot ?? {
              effective_peaks: [],
              confirmed_index: null,
              analysis_inputs_hash: "",
            });
        const out: ComparisonMemberInput = {
          exposure_id: m.exposure_id,
          display_order: m.display_order,
          band_height: m.band_height,
          y_offset: m.y_offset,
          normalization: m.normalization,
          snapshot,
        };
        if (m.id !== undefined) out.id = m.id;
        if (m.color_override !== undefined) out.color_override = m.color_override;
        if (m.label_override !== undefined) out.label_override = m.label_override;
        if (m.q_window_min !== undefined) out.q_window_min = m.q_window_min;
        if (m.q_window_max !== undefined) out.q_window_max = m.q_window_max;
        if (m.peak_display !== undefined) out.peak_display = m.peak_display;
        return out;
      });
      // useSaveComparison flat-spreads the input into the SaveComparisonBody;
      // see saveComparison mutator's `request: (p) => api.saveComparison(...)`.
      const payload: SaveComparisonBody & { id?: number } = {
        title: d.title,
        members,
      };
      if (d.id !== undefined) payload.id = d.id;
      if (d.description !== "") payload.description = d.description;
      if (d.baseHash !== undefined) payload.expected_content_hash = d.baseHash;
      // Phase 11 — fork lineage rides through to POST /api/comparisons. Both
      // fields ride together (or not at all) per backend contract; the UI
      // factory `fromComparisonAsFork` always sets both when populating a fork.
      if (d.forkedFromId !== undefined) payload.forked_from_id = d.forkedFromId;
      if (d.forkedAtHash !== undefined) payload.forked_at_hash = d.forkedAtHash;
      // C-4 — forward author's view choices so the server persists them.
      // undefined (never set) ⇒ sent as null (backend clears / uses default);
      // a value ⇒ sent as-is (backend stores it).
      payload.view_grouping_mode    = d.viewGroupingMode    ?? null;
      payload.view_show_peak_ticks  = d.viewShowPeakTicks   ?? null;
      payload.view_show_peak_labels = d.viewShowPeakLabels  ?? null;
      return payload;
    },
    [qc],
  );

  const setDraftForkOf = useAppState((s) => s.setDraftForkOf);

  const handleSave = useCallback(async () => {
    if (draft === null) return;
    if (draft.members.length === 0) return;
    // Guard against duplicate triggers during the async prefetch window.
    // save.isPending is false until save.mutate() fires, so a second
    // Cmd+Enter while awaiting would start a parallel round. Set the
    // in-flight ref early so the keyboard handler and button both reject.
    if (pendingSubmitRef.current) return;

    // ── Compare UX C-14 — editing-as-fork-of save path ──────────────────
    // A NON-author saving a draft on someone else's comparison must NOT
    // overwrite the original. Prompt for a fork title, morph the draft into
    // a fork (clear id/baseHash, set lineage), then submit via the create
    // path. The morphed draft's `id === undefined` routes the mutator to
    // POST /api/comparisons with no `expected_content_hash`.
    if (compareMode.kind === "editing-as-fork-of") {
      const baseTitle = comparisonQ.data?.title ?? "comparison";
      // window.prompt note: Playwright auto-dismisses dialogs by default.
      // Any e2e covering this flow MUST register a dialog handler before
      // the gesture: page.on("dialog", d => d.accept("My fork")).
      const proposed = window.prompt("Title for your fork:", `Copy of ${baseTitle}`);
      if (proposed === null) return; // cancelled — leave the draft intact
      setDraftForkOf({
        newTitle: proposed.trim() === "" ? `Copy of ${baseTitle}` : proposed,
        sourceId: compareMode.parentId,
        sourceHash: draft.baseHash ?? "",
      });
      pendingSubmitRef.current = true;
      // Re-read the MORPHED draft from the store before building the
      // payload. Reusing a payload built from the pre-morph `draft` closure
      // would carry the original id and the create-path route would never
      // fire (and the two back-to-back ops would collide on one
      // `idempotent_responses` row). Build fresh, post-morph.
      const fresh = useAppState.getState().activeDraft;
      const exposureIds = (fresh?.members ?? [])
        .map((m) => m.exposure_id)
        .filter((id): id is number => id !== null);
      try {
        await prefetchColdMembers(exposureIds, qc);
        const payload = buildSavePayload(fresh);
        if (payload === null) {
          pendingSubmitRef.current = false;
          return;
        }
        save.mutate(payload);
      } catch {
        pendingSubmitRef.current = false;
      }
      return;
    }

    pendingSubmitRef.current = true;

    // Warm the four cache keys computeMemberSnapshot reads (#49) before
    // computing snapshots — without this, never-visited members land with
    // analysis_inputs_hash = "" and the server marks them stale on the
    // next view fold. Per-key cold detection (#93): each key is checked
    // independently so an exposure with three warm keys and one cold key
    // only refetches the missing one.
    const exposureIds = draft.members
      .map((m) => m.exposure_id)
      .filter((id): id is number => id !== null);

    try {
      await prefetchColdMembers(exposureIds, qc);
      const payload = buildSavePayload(draft);
      if (payload === null) {
        pendingSubmitRef.current = false;
        return;
      }
      save.mutate(payload);
    } catch {
      // Prefetch or mutate failed — release the guard so the user can retry.
      pendingSubmitRef.current = false;
    }
  }, [draft, qc, save, compareMode, comparisonQ.data, setDraftForkOf, buildSavePayload]);

  // Post-success navigation. Reading `save.data` (the response) lets us
  // navigate to /experiments/:eid/compare/:newId for create flow. We gate
  // on `pendingSubmitRef` so only a save we just initiated triggers nav.
  useEffect(() => {
    if (!save.isSuccess) return;
    if (!pendingSubmitRef.current) return;
    pendingSubmitRef.current = false;
    const response = save.data as Comparison | undefined;
    discardDraft();
    if (response && typeof response.id === "number") {
      goToReview(response.id);
    } else {
      goToList();
    }
  }, [save.isSuccess, save.data, discardDraft, goToReview, goToList]);

  // Release the guard on mutation error so the user can retry.
  // save.mutate() is fire-and-forget; errors surface via save.error.
  useEffect(() => {
    if (save.error) pendingSubmitRef.current = false;
  }, [save.error]);

  // Phase 13 Task 13.4 — keyboard shortcuts:
  //   Esc            → cancel (return to review or list)
  //   Cmd/Ctrl+Enter → submit (Save), if save isn't already disabled
  // Listener is attached while the edit body is mounted so the shortcut only
  // fires during an edit session. The save handler reads `draft` through the
  // closure already; we re-check `members.length === 0` to mirror the
  // button's disabled state and avoid sending an empty payload.
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        // Don't intercept Esc when a modal (picker, conflict, nav) is open
        // — those have their own Esc handling. `:not([hidden])` filters out
        // dialogs that exist in the DOM but are not currently presented (a
        // future regression where a hidden modal silently swallows Esc).
        const openDialog = document.querySelector('[role="dialog"]:not([hidden])');
        if (openDialog) return;
        e.preventDefault();
        handleCancel();
        return;
      }
      if ((e.metaKey || e.ctrlKey) && e.key === "Enter") {
        if (draft === null || draft.members.length === 0) return;
        if (save.isPending) return;
        e.preventDefault();
        handleSave();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [handleCancel, handleSave, draft, save.isPending]);

  // ── Compare UX C-13 — edit-mode header wiring ────────────────────────────
  // The edit header uses the shared CompareTitleStrip / CompareStatusSurface
  // / CompareToolbar components plus a SavePill (issue #139), mirroring the
  // C-12 review-mode wiring. `currentUserId` / `compareMode` are declared
  // above `handleSave` (C-14).
  const forksQ = useComparisonForks(id);
  // Pre-resolved fork links for the toolbar's ⋯-more dropdown (C-17).
  const forksList = useMemo(
    () => (forksQ.data ?? []).map((f) => ({
      id: f.id,
      title: f.title,
      href: comparePath({ scope, eid, id: f.id }),
    })),
    [forksQ.data, scope, eid],
  );
  const deleteMut = useDeleteComparison();

  // Toolbar overflow-menu handlers. Edit mode has no live comparison id for
  // the create flow, so Copy-link / Fork / Delete only do meaningful work
  // when editing an existing comparison.
  const handleCopyLink = useCallback(() => {
    void navigator.clipboard?.writeText(window.location.href);
  }, []);

  const handleDelete = useCallback(() => {
    if (id === undefined) return;
    deleteMut.mutate({ id });
    discardDraft();
    goToList();
  }, [deleteMut, id, discardDraft, goToList]);

  // Edit mode is itself the fork-authoring surface; "Fork" from here just
  // returns to the create flow with the current draft intact.
  const handleFork = useCallback(() => {
    navigate(comparePath({ scope, eid, isNew: true }));
  }, [navigate, scope, eid]);

  // Grouping-mode write — mirrors the review body; auto-creates an empty
  // draft when none is active (spec §6.4 viewer escape hatch).
  const setDraftViewGroupingMode = useAppState((s) => s.setDraftViewGroupingMode);

  // SavePill renders only when the draft is dirty + non-empty; the legacy
  // Save button was disabled (not hidden) for an empty draft. Treat a draft
  // with at least one member as "dirty" so the pill surfaces the save
  // affordance, keeping C-13 behaviour-equivalent to the old button.
  const saveDirty = (draft?.members.length ?? 0) > 0;

  // ── Compare UX C-15 Step 4 — empty-state handlers ─────────────────────────
  // `handleAddTraces` MUST branch on `mode.kind`. In edit mode the draft is
  // already a draft surface, so the gesture just focuses the picker's search
  // input — the interim picker chip. In a hypothetical viewing branch it
  // would seed a draft from the current comparison via `fromComparison`
  // (the real factory; `loadDraftFromComparison` is the matching action) so
  // the next save updates the comparison rather than minting a blank one.
  // The edit body only ever renders while a draft is active or on `/new`, so
  // the viewing branch is defensive — but kept faithful to the C-15 spec.
  const handleAddTraces = useCallback(() => {
    if (compareMode.kind === "viewing" || compareMode.kind === "viewing-stale") {
      if (comparisonQ.data !== undefined) {
        loadDraftFromComp(comparisonQ.data, qc);
      }
      return;
    }
    // creating-* / editing-* : already a draft surface — open the interim
    // picker chip by focusing its search input.
    pickerSearchRef.current?.focus();
  }, [compareMode.kind, comparisonQ.data, loadDraftFromComp, qc]);

  // Drop-target wiring for external exposure drags (§7.1). External-add is
  // fully Phase 2; for now the drop handler is a thin route to the interim
  // picker — `preventDefault` on dragover marks the region as a valid drop
  // target, and a drop focuses the picker search.
  const handleExposureDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.dataTransfer.dropEffect = "copy";
  }, []);
  const handleExposureDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    pickerSearchRef.current?.focus();
  }, []);

  // Phase 6 — wire MultiTracePlot. Members reflect the live draft (Zustand);
  // traces fetched in parallel via `useMemberTraces`. The q-axis zoom domain
  // is keyed per comparison id so toggling between review/edit for the SAME
  // comparison preserves the zoom but switching to a DIFFERENT comparison
  // doesn't bleed (zoom ranges can differ by orders of magnitude).
  // Unsaved drafts (create flow, no id yet) use the `0` sentinel — autoincrement
  // comparison ids start at 1 so collision is impossible.
  const xDomainKey = id ?? 0;
  const xDomain = useAppState((s) => s.compareXDomains[xDomainKey] ?? null);
  const setCompareXDomain = useAppState((s) => s.setCompareXDomain);
  const setXDomain = useCallback(
    (d: [number, number] | null) => setCompareXDomain(xDomainKey, d),
    [setCompareXDomain, xDomainKey],
  );
  const plotMembers = useMemo<SeriesMember[]>(
    () => (draft?.members ?? []).map(draftToMember),
    [draft?.members],
  );

  // Phase 8.1 — peak click cycle. Maps a `MultiTracePlot` callback
  // (memberId, peakId, altKey) → Zustand cycle action by looking up the
  // member index from its synthetic id (draftToMember mints stable ids).
  const cyclePeak = useAppState((s) => s.cyclePeakDisplayForMember);
  const handlePeakClick = useCallback(
    (memberId: number, peakId: number, altKey: boolean) => {
      const idx = plotMembers.findIndex((m) => m.id === memberId);
      if (idx < 0) return;
      cyclePeak(idx, peakId, altKey);
    },
    [plotMembers, cyclePeak],
  );

  // Materialize the per-member peak_display map from draft state so the
  // plot's optimistic state survives between server-side updates.
  const peakDisplayByMemberId = useMemo(() => {
    const m = new Map<number, { hidden: number[]; labeled: number[] }>();
    for (const dm of draft?.members ?? []) {
      const memberId = dm.id ?? -(dm.display_order + 1);
      if (dm.peak_display) m.set(memberId, dm.peak_display);
    }
    return m;
  }, [draft?.members]);
  const exposureIds = useMemo(
    () => plotMembers.flatMap((m) => (m.exposure_id !== null ? [m.exposure_id] : [])),
    [plotMembers],
  );
  const traces = useMemberTraces(exposureIds);
  const tracesLoading = useMemberTracesLoading(exposureIds);
  const resetBandHeights = useAppState((s) => s.resetBandHeights);

  // Issue #69 — hydrate per-member exposure + sample rows so the gutter
  // resolves human-readable labels and `sampleIdFor` resolves grouping
  // colors. Picker pre-warming (#49) covers the picker-add path; this
  // covers the deep-link cold load (`/experiments/:eid/compare/:id`),
  // where the comparison fetch hydrates the draft but never the
  // per-exposure rows. Mirror of the review body's #61 + #52 fix.
  const exposures = useMemberExposures(exposureIds);
  const sampleIds = useMemo(() => {
    const ids = new Set<number>();
    for (const e of exposures.values()) ids.add(e.sample_id);
    return Array.from(ids).sort((a, b) => a - b);
  }, [exposures]);
  const samples = useMemberSamples(sampleIds);

  // Boneyard fixture mirrors the dual-column plot+gutter geometry.
  const editPlotFixture = useMemo(
    () => (
      <div className="flex flex-row flex-1 min-h-0 gap-3">
        <div className="flex-1 min-w-0 border border-border/40 rounded h-full" />
        <div className="w-[320px] shrink-0 border border-border/40 rounded h-full" />
      </div>
    ),
    [],
  );

  // Phase 9 gap-fix — line-stroke coloring grouping mode + sample-id resolver
  // for edit mode. Mirrors the review body so toggling the grouping mode in
  // sister tabs / re-entering review keeps trace colors consistent. Reads
  // from the per-member exposures Map (hydrated by useMemberExposures
  // above) instead of `qc.getQueryData` so the resolver re-evaluates when
  // the cache settles — see issue #69.
  //
  // C-4 Step 0: groupingMode is now resolved via effectiveGroupingMode so the
  // draft's viewGroupingMode takes precedence over the server record.
  const groupingMode = effectiveGroupingMode(draft, comparisonQ.data);
  const sampleIdFor = useCallback(
    (m: SeriesMember): number | null => {
      if (m.exposure_id === null) return null;
      return exposures.get(m.exposure_id)?.sample_id ?? null;
    },
    [exposures],
  );

  // Issue #69 — same fallback chain as the review body. Both share
  // `resolveDisplayLabels` so the chain stays in lockstep.
  const displayLabelByMemberId = useMemo(
    () => resolveDisplayLabels(plotMembers, exposures, samples),
    [plotMembers, exposures, samples],
  );

  // Track the plot column's height so the edit-mode gutter aligns pixel-for-
  // pixel with the plot bands (both consumers feed `computeYBands` the same
  // ratios + height).
  //
  // Re-attach the observer on `tracesLoading` flips: while the Skeleton
  // fallback is showing, `plotColRef.current` is null and the effect bails
  // at the early-out. Once Skeleton swaps in the real children, the ref
  // attaches and we need to fire again — same shape as the review fix
  // (issue #51, see PR #59 [plotLoading]). `plotMembers.length` alone
  // doesn't gate on Skeleton, so the gutter stayed at height 0 until any
  // other re-render happened to coincide with the ref attaching.
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
  }, [tracesLoading]);

  const members = draft?.members ?? [];

  // Edit-mode center pane — title strip + description + plot host. Lives
  // in the WorkspaceGrid center slot; the right slot holds the interim
  // ComparisonPickerPanel chip.
  const editCenter = (
    <div className="flex-1 min-h-0 flex flex-col p-4 gap-3">
      <div data-testid="compare-edit-header" className="flex flex-col gap-2">
        <CompareTitleStrip
          title={draft?.title ?? ""}
          // Editable (non-readonly): `""` keeps the description row visible
          // with the "Add a description…" placeholder so it can be edited.
          description={draft?.description ?? ""}
          memberCount={members.length}
          // Editing your own draft — the byline reads "by you".
          authorUsername={null}
          isCurrentUserAuthor
          lastEventAt={comparisonQ.data?.last_event_at ?? null}
          forkedFromTitle={comparisonQ.data?.forked_from_title ?? null}
          forkedFromHref={null}
          onTitleChange={setDraftTitle}
          onDescChange={setDraftDescription}
        />
        <CompareStatusSurface
          needsReview={null}
          serverUpdate={null}
          savedAt={null}
        />
        <CompareToolbar
          groupingControl={
            <GroupingModeToggle
              mode={groupingMode}
              onChange={setDraftViewGroupingMode}
            />
          }
          // Edit mode has never carried peak-annotation toggles (the edit
          // plot doesn't read `showPeakTicks`/`showPeakLabels`); the
          // annotation slot instead hosts the edit-only "Reset heights"
          // control, which has no home in the shared component set.
          annotationControl={
            <button
              type="button"
              data-testid="compare-edit-reset-heights"
              onClick={resetBandHeights}
              disabled={members.length === 0}
              className="px-2 py-0.5 rounded border border-border text-fg-muted text-xs
                         disabled:opacity-50"
              title="Reset all band heights to default"
            >
              Reset heights
            </button>
          }
          forksList={forksList}
          onCopyLink={handleCopyLink}
          onDelete={handleDelete}
          onDiscardChanges={handleDiscard}
          onFork={handleFork}
          exportControl={null}
          saveControl={
            <span className="inline-flex items-center gap-2">
              <button
                type="button"
                data-testid="compare-edit-cancel"
                onClick={handleCancel}
                title="Cancel (Esc)"
                className="px-3 py-1 rounded border border-border text-sm"
              >
                Cancel
              </button>
              <SavePill
                dirty={saveDirty}
                mode={compareMode}
                onSave={handleSave}
                isSaving={save.isPending}
              />
            </span>
          }
        />
      </div>
      <div
        data-testid="compare-edit-plot-host"
        className="flex-1 min-h-0 flex flex-col gap-2"
      >
        {plotMembers.length === 0 ? (
          // ── Compare UX C-15 Step 4 — §7.1 empty state ──────────────────────
          // Three-element empty block: headline, "+ Add traces" button, and
          // a drag hint. The region is also a drop target for external
          // exposure drags. Carries `compare-empty-state` (the canonical
          // C-15 testid) plus `compare-edit-plot-empty` continuity is via the
          // inner button keeping `compare-edit-add-traces`.
          <div
            data-testid="compare-empty-state"
            onDragOver={handleExposureDragOver}
            onDrop={handleExposureDrop}
            className="flex-1 flex flex-col items-center justify-center
                       border border-border/40 rounded text-fg-muted text-sm gap-3"
          >
            <h2 className="text-fg font-medium">No traces yet.</h2>
            <button
              type="button"
              data-testid="compare-edit-add-traces"
              data-interactable="button"
              onClick={handleAddTraces}
              className="px-3 py-1 rounded border border-border text-fg text-sm
                         hover:bg-bg-elevated"
            >
              + Add traces
            </button>
            <p>Or drag exposures from the picker.</p>
          </div>
        ) : (
          <Skeleton
            name="compare-edit-plot"
            className="flex flex-row flex-1 min-h-0 gap-3"
            loading={tracesLoading}
            stagger={50}
            transition={200}
            fixture={editPlotFixture}
            fallback={<div className="flex-1 flex items-center justify-center"><HintText>Loading traces…</HintText></div>}
          >
            {/*
              Compare UX E-4 — the inter-row resize gaps in `MemberMetaGutter`
              publish the band-above id; the per-band overlays in
              `MultiTracePlot` subscribe to tint accent on hover/drag. The
              two are otherwise unrelated siblings, so the coupling rides a
              minimal context scoped to this edit-mode plot+gutter pair.
            */}
            <ActiveBandProvider>
              <div ref={plotColRef} className="flex-1 min-w-0">
                <MultiTracePlot
                  members={plotMembers}
                  traces={traces}
                  xDomain={xDomain}
                  onXDomain={setXDomain}
                  peakDisplayByMemberId={peakDisplayByMemberId}
                  onPeakClick={handlePeakClick}
                  groupingMode={groupingMode}
                  sampleIdFor={sampleIdFor}
                />
              </div>
              <div
                className="w-[320px] shrink-0 relative"
                data-testid="compare-edit-gutter"
              >
                <MemberMetaGutter
                  members={plotMembers}
                  panelHeight={panelHeight}
                  mode="edit"
                  displayLabelByMemberId={displayLabelByMemberId}
                />
              </div>
            </ActiveBandProvider>
          </Skeleton>
        )}
      </div>
    </div>
  );

  return (
    <div
      data-testid="compare-page-edit"
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
        center={editCenter}
        right={
          /* Interim picker chip — Phase 2 will move this to a right-rail tab. */
          <ComparisonPickerPanel experimentId={eid} searchInputRef={pickerSearchRef} />
        }
        slotClassName={{
          // ComparisonSidebar uses `flex-1`, so the slot needs `display:flex`.
          left:   "flex flex-col min-h-[400px]",
          center: "flex flex-col min-h-[640px]",
          right:  "flex flex-col min-h-[400px]",
        }}
      />
    </div>
  );
}
