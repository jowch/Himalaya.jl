/**
 * ComparePageEdit — edit-mode shell with save flow (Plan §Phase 4, Task 4.4).
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare/new        — empty draft (create flow)
 *   /experiments/:eid/compare/:id/edit   — edit existing comparison
 *
 * Owns:
 *   - Title + description inputs (controlled, bound to Zustand draft).
 *   - Save / Cancel / Discard buttons.
 *   - Snapshot recompute at submit time per Plan §Task 4.3.
 *
 * The plot, member panel, drag handles, picker chip, and conflict modal
 * land in Phases 5–12; this file only mounts placeholders for them.
 *
 * Save flow (Plan §Task 4.4):
 *   - Disabled when `draft.members.length === 0`.
 *   - Computes a fresh `MemberSnapshot` per member at submit (cache-derived).
 *   - Calls `useSaveComparison`. The mutator routes to POST /api/comparisons
 *     (create) or POST /api/comparisons/:id/submit (update) based on the
 *     `id` field; `expected_content_hash` rides from `draft.baseHash`.
 *   - On success: clear draft, navigate to /experiments/:eid/compare/:newId.
 *   - On 409 ConflictError: leave the draft intact. The conflict modal in
 *     Phase 12 picks up the typed throw via `save.error`.
 *
 * The hydration effect aligns the in-memory draft with the URL on every
 * mount: `/new` ⇒ start a fresh empty draft if one isn't already active;
 * `/:id/edit` ⇒ load from the comparison fetch (with cache-derived snapshot
 * recovery) if the draft isn't already this id.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { Skeleton } from "boneyard-js/react";
import { ComparisonSidebar } from "../components/ComparisonSidebar";
import { ComparisonPicker } from "../components/ComparisonPicker";
import { MultiTracePlot } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { HintText } from "../components/ui";
import { useAppState } from "../state";
import {
  useSaveComparison, useComparison, useMemberTraces, useMemberTracesLoading, queryKeys,
} from "../queries";
import { computeMemberSnapshot } from "../lib/comparison/snapshot";
import type {
  Comparison, ComparisonMember, ComparisonMemberInput, Exposure, SaveComparisonBody,
} from "../api";
import type { DraftMember } from "../lib/comparison/draft";

/**
 * Convert a draft member into a ComparisonMember-shaped object suitable for
 * `MultiTracePlot`. Unsaved drafts have `id = undefined`; we substitute a
 * stable negative synthetic id keyed by display_order so the plot's per-member
 * keying stays consistent across re-renders. Snapshot can also be undefined
 * mid-edit; the plot tolerates a null snapshot (no peaks rendered).
 */
function draftToMember(d: DraftMember): ComparisonMember {
  return {
    id: d.id ?? -(d.display_order + 1),
    comparison_id: 0,
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

export function ComparePageEdit(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const navigate = useNavigate();
  const qc = useQueryClient();

  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;

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
  const comparisonQ = useComparison(id);
  useEffect(() => {
    if (id === undefined) {
      // /new: start an empty draft if there isn't one (or one that's tied
      // to a different comparison id). Don't clobber a fresh draft mid-edit.
      if (draft === null || draft.id !== undefined) {
        startNewDraft();
      }
    } else if (comparisonQ.data && (draft === null || draft.id !== id)) {
      // /:id/edit: load from the server fetch with cache-derived snapshots.
      loadDraftFromComp(comparisonQ.data, qc);
    }
    // Intentionally exclude `draft` from deps — we only re-evaluate when
    // the URL id or the fetched comparison changes. Reading `draft` inside
    // the effect for the comparison-id check is fine; including it would
    // re-run the effect after the load and re-load again.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [id, comparisonQ.data, startNewDraft, loadDraftFromComp, qc]);

  const save = useSaveComparison();
  const pendingSubmitRef = useRef(false);

  const goToReview = useCallback(
    (newId: number) => {
      if (eid !== undefined) navigate(`/experiments/${eid}/compare/${newId}`);
      else navigate(`/compare/all`);
    },
    [navigate, eid],
  );

  const goToList = useCallback(() => {
    if (eid !== undefined) navigate(`/experiments/${eid}/compare`);
    else navigate(`/compare/all`);
  }, [navigate, eid]);

  const handleCancel = useCallback(() => {
    if (id !== undefined) {
      // Cancel from edit-existing → return to that comparison's review page.
      if (eid !== undefined) navigate(`/experiments/${eid}/compare/${id}`);
      else navigate(`/compare/all`);
    } else {
      // Cancel from create → return to list.
      goToList();
    }
  }, [navigate, id, eid, goToList]);

  const handleDiscard = useCallback(() => {
    discardDraft();
    goToList();
  }, [discardDraft, goToList]);

  const handleSave = useCallback(() => {
    if (draft === null) return;
    if (draft.members.length === 0) return;
    // Compute a fresh snapshot per member at submit time (Plan §Task 4.3).
    const members: ComparisonMemberInput[] = draft.members.map((m) => {
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
      title: draft.title,
      members,
    };
    if (draft.id !== undefined) payload.id = draft.id;
    if (draft.description !== "") payload.description = draft.description;
    if (draft.baseHash !== undefined) payload.expected_content_hash = draft.baseHash;
    // Phase 11 — fork lineage rides through to POST /api/comparisons. Both
    // fields ride together (or not at all) per backend contract; the UI
    // factory `fromComparisonAsFork` always sets both when populating a fork.
    if (draft.forkedFromId !== undefined) payload.forked_from_id = draft.forkedFromId;
    if (draft.forkedAtHash !== undefined) payload.forked_at_hash = draft.forkedAtHash;
    // Mark this submit as "in-flight" so the post-success effect knows to
    // navigate. Without the ref, an already-saved success state on remount
    // would re-fire navigation.
    pendingSubmitRef.current = true;
    save.mutate(payload);
  }, [draft, qc, save]);

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

  // Phase 13 Task 13.4 — keyboard shortcuts:
  //   Esc            → cancel (return to review or list)
  //   Cmd/Ctrl+Enter → submit (Save), if save isn't already disabled
  // Listener is attached to the page-edit shell so the shortcut only fires
  // while the user is on the edit page. The save handler reads `draft`
  // through the closure already; we re-check `members.length === 0` to
  // mirror the button's disabled state and avoid sending an empty payload.
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        // Don't intercept Esc when a modal (picker, conflict, nav) is open
        // — those have their own Esc handling. The presence of an open
        // dialog with role="dialog" is the cleanest cross-modal probe.
        const openDialog = document.querySelector('[role="dialog"]');
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

  const scope: "all" | "experiment" = eid !== undefined ? "experiment" : "all";

  // Phase 5 Task 5.2 — local state for the picker open/close. Picker open
  // state is purely client-side, no need for Zustand.
  const [pickerOpen, setPickerOpen] = useState(false);

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
  const plotMembers = useMemo<ComparisonMember[]>(
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
  // for edit mode. Mirrors ComparePage so toggling the grouping mode in
  // sister tabs / re-entering review keeps trace colors consistent. Reads
  // the TanStack `exposure` cache directly; cache misses → ORPHAN_FALLBACK
  // until the fetch settles.
  const groupingMode = useAppState((s) => s.groupingMode);
  const sampleIdFor = useCallback(
    (m: ComparisonMember): number | null => {
      if (m.exposure_id === null) return null;
      const exposure = qc.getQueryData<Exposure>(queryKeys.exposure(m.exposure_id));
      return exposure?.sample_id ?? null;
    },
    [qc],
  );

  // Track the plot column's height so the edit-mode gutter aligns pixel-for-
  // pixel with the plot bands (both consumers feed `computeYBands` the same
  // ratios + height).
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
  }, [plotMembers.length]);
  // Phase 11 wires Edit/Fork visibility against current_user vs. created_by;
  // for now we surface the testid so downstream tests can target it once
  // the gating exists. The button is hidden from the rendered tree until
  // then to keep the edit shell clean.

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
      <section className="card overflow-hidden flex-1 min-h-0 flex flex-col p-4 gap-3">
        <div className="flex items-center gap-2">
          <input
            data-testid="compare-edit-title"
            type="text"
            placeholder="Comparison title"
            value={draft?.title ?? ""}
            onChange={(e) => setDraftTitle(e.target.value)}
            className="flex-1 bg-transparent border border-border rounded px-2 py-1 text-base"
          />
          <button
            type="button"
            data-testid="comparison-save"
            onClick={handleSave}
            disabled={(draft?.members.length ?? 0) === 0 || save.isPending}
            title="Save (Cmd+Enter)"
            className="px-3 py-1 rounded bg-accent text-bg disabled:opacity-50 text-sm font-medium"
          >
            {save.isPending ? "Saving…" : "Save"}
          </button>
          <button
            type="button"
            data-testid="comparison-cancel"
            onClick={handleCancel}
            title="Cancel (Esc)"
            className="px-3 py-1 rounded border border-border text-sm"
          >
            Cancel
          </button>
          <button
            type="button"
            data-testid="comparison-discard"
            onClick={handleDiscard}
            className="px-3 py-1 rounded border border-border text-fg-muted text-sm"
          >
            Discard draft
          </button>
          <button
            type="button"
            data-testid="compare-edit-reset-heights"
            onClick={resetBandHeights}
            disabled={(draft?.members.length ?? 0) === 0}
            className="px-3 py-1 rounded border border-border text-fg-muted text-sm
                       disabled:opacity-50"
            title="Reset all band heights to default"
          >
            Reset heights
          </button>
        </div>
        <textarea
          data-testid="compare-edit-description"
          placeholder="Description (optional)"
          value={draft?.description ?? ""}
          onChange={(e) => setDraftDescription(e.target.value)}
          className="bg-transparent border border-border rounded px-2 py-1 text-sm resize-none h-16"
        />
        <div
          data-testid="compare-edit-plot-host"
          className="flex-1 min-h-0 flex flex-col gap-2"
        >
          {plotMembers.length === 0 ? (
            <div
              data-testid="compare-edit-plot-empty"
              className="flex-1 flex flex-col items-center justify-center
                         border border-border/40 rounded text-fg-muted text-sm gap-3"
            >
              <div>No traces yet — add some to get started.</div>
              <button
                type="button"
                data-testid="compare-edit-add-traces"
                onClick={() => setPickerOpen(true)}
                className="px-3 py-1 rounded border border-border text-fg text-sm
                           hover:bg-bg-elevated"
              >
                + Add traces
              </button>
            </div>
          ) : (
            <>
              <Skeleton
                name="compare-edit-plot"
                className="flex flex-row flex-1 min-h-0 gap-3"
                loading={tracesLoading}
                stagger={50}
                transition={200}
                fixture={editPlotFixture}
                fallback={<div className="flex-1 flex items-center justify-center"><HintText>Loading traces…</HintText></div>}
              >
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
                  />
                </div>
              </Skeleton>
              <div className="flex justify-end">
                <button
                  type="button"
                  data-testid="compare-edit-add-traces"
                  onClick={() => setPickerOpen(true)}
                  className="px-3 py-1 rounded border border-border text-fg text-sm
                             hover:bg-bg-elevated"
                >
                  + Add traces
                </button>
              </div>
            </>
          )}
        </div>
      </section>

      <ComparisonPicker
        isOpen={pickerOpen}
        onClose={() => setPickerOpen(false)}
        experimentId={eid}
      />
    </div>
  );
}
