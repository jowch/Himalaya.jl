import { create } from "zustand";
import { persist } from "zustand/middleware";
import type { QueryClient } from "@tanstack/react-query";
import type { Comparison, ConflictError } from "./api";
import {
  emptyDraft,
  loadDraftFromSession,
  persistDraftToSession,
  type ActiveDraft,
  type ActiveDraftSlot,
  type DraftMember,
} from "./lib/comparison/draft";
import {
  fromComparison,
  fromComparisonAsFork,
  memberFromNewExposure,
} from "./lib/comparison/draftFactories";
import { cyclePeakDisplay } from "./lib/comparison/peakCycle";
import { emitReplaceNext } from "./lib/url/emitMode";

export const LS_KEY = "himalaya-ui:state";

export type PageId = "index" | "compare" | "inspect";
export type ThemeId = "dark" | "light";
export type NavModalStep = "experiment" | "sample";

export type StaleUrlContext =
  | {
      kind: "not_found";
      missing: "experiment" | "sample" | "exposure";
      missing_value: string;
      experiment_resolved: { id: number; name: string } | undefined;
      sample_resolved: { id: number; name: string } | undefined;
    }
  | { kind: "unknown_path"; raw: string };

export type RecoverOpts = {
  step: NavModalStep;
  experimentId: number | undefined;
  sampleId: number | undefined;
  openModal?: boolean; // default true; row "exposure" passes false
};

export type ResolveSuccessSlots = {
  page: PageId;
  experimentId: number | undefined;
  sampleId: number | undefined;
  exposureId: number | undefined;
};

export interface AppState {
  // persisted
  username: string | undefined;
  firstName: string | undefined;
  lastName: string | undefined;
  activeExperimentId: number | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  activePage: PageId;
  tutorialSeen: boolean;
  theme: ThemeId;

  // ephemeral (not persisted to localStorage)
  hoveredIndexId: number | undefined;
  hoveredPeakId: number | undefined;
  navModalOpen: boolean;
  navModalStep: NavModalStep;
  // Speculative builder: when non-null, the modal is open for this exposure.
  // All builder form state (phase, anchor peak, ratio) is local to the
  // SpeculativeBuilder component — only the open/close gate lives in store
  // because PhasePanel needs to mount/unmount the modal.
  speculativeBuilder: { exposureId: number } | null;
  /**
   * Compare-page q-axis zoom domains, keyed per comparison id. Per-tab UI
   * state — not persisted. Missing entry / `null` value = full data range.
   *
   * Keying per-comparison preserves zoom across review/edit toggles for
   * the SAME comparison while isolating different comparisons (different
   * comparisons can have q-ranges differing by orders of magnitude — a
   * shared single slot caused B to inherit A's zoom and look broken).
   *
   * The unsaved-draft case (create flow, no id yet) uses key `0` as a
   * sentinel — autoincrement comparison ids start at 1, so collision is
   * impossible.
   */
  compareXDomains: Record<number, [number, number] | null>;

  /**
   * Compare-page draft slot (Plan §Phase 4, Task 4.3). Single slot — only
   * one comparison can be in edit mode at a time per tab. Mirrored to
   * sessionStorage with a schema version (see `lib/comparison/draft.ts`).
   * Tab close loses the draft, which is acceptable for v1 per the spec.
   */
  activeDraft: ActiveDraftSlot;

  /**
   * Compare-page review-mode annotation toggles (Plan §Phase 9, Task 9.3).
   * Both default to `true`. Per-tab viewing preference — neither persisted
   * on the comparison nor in storage. Hidden in edit mode (everything
   * renders so the user can edit).
   */
  showPeakTicks: boolean;
  showPeakLabels: boolean;

  /**
   * Compare-page hover/click-to-pin highlight target (Plan §Phase 9,
   * Task 9.5). When set, `MemberTraceLayer` recolors that member's
   * snapshotted index peaks to the phase color; non-index peaks stay black.
   * Mirrors the `hoveredIndexId` single-setter pattern from the Index page.
   * Cleared on page navigation, edit-mode entry, and member removal.
   */
  highlightedCompareMemberId: number | undefined;

  /**
   * Compare-page conflict modal slot (Plan §Phase 12). When non-null, the
   * `ConflictModal` mounted at `App.tsx` opens, rendering the server's
   * `current_state` (frozen at conflict time) side-by-side with the local
   * draft. Set by `useSaveComparison` whenever the typed `ConflictError`
   * surfaces; cleared by Discard / Overwrite-success / Fork / Esc.
   *
   * NOT persisted — a 409 is a tab-local UX concern, and replaying it
   * across reloads would resurrect a stale conflict whose underlying
   * server state has likely moved on.
   *
   * Re-callability invariant: setting a fresh `ConflictError` while the
   * modal is open replaces the slot rather than stacking. This is the
   * second-409 race path — the modal stays mounted but its rendered
   * server-state panel updates to the new `current_state`.
   */
  pendingConflict: ConflictError | null;

  /**
   * Permalink URL handling slots (spec §4.4 + §6).
   * Both ephemeral — not persisted. `staleUrlContext` is non-null when the
   * current URL points to a slug that doesn't resolve (404 from
   * `/api/resolve` or unknown path). `resolving` is true while the
   * URL→state resolve fetch is in flight.
   */
  staleUrlContext: StaleUrlContext | null;
  resolving: boolean;

  // setters
  setUsername: (name: string) => void;
  setUser: (u: { username: string; firstName?: string | undefined; lastName?: string | undefined }) => void;
  setActiveExperiment: (id: number | undefined) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
  setHoveredPeak: (id: number | undefined) => void;
  setActivePage: (page: PageId) => void;
  setTutorialSeen: (seen: boolean) => void;
  setTheme: (theme: ThemeId) => void;
  openNavModal: (step?: NavModalStep) => void;
  closeNavModal: () => void;
  setNavModalStep: (step: NavModalStep) => void;
  clearUsername: () => void;
  openSpeculativeBuilder: (exposureId: number) => void;
  closeSpeculativeBuilder: () => void;
  setCompareXDomain: (id: number, d: [number, number] | null) => void;

  // Compare-draft actions
  startNewDraft: () => void;
  /**
   * Start a fork-flavored draft pre-populated from a parent comparison
   * (Plan §Phase 11, Task 11.2). Members come from the parent (with ids
   * dropped so they INSERT under the new comparison) and the parent's
   * lineage rides on the draft so the eventual `POST /api/comparisons`
   * carries `forked_from_id` + `forked_at_hash`.
   */
  startForkDraft: (comparison: Comparison, qc: QueryClient) => void;
  loadDraftFromComparison: (comparison: Comparison, qc: QueryClient) => void;
  setDraftTitle: (title: string) => void;
  setDraftDescription: (description: string) => void;
  /**
   * Morph the active draft into a fork (Compare UX C-14). Called when a
   * NON-author saves a draft on someone else's comparison: clears `id` +
   * `baseHash` (so the next submit routes to the create path) and records
   * the parent lineage (`forkedFromId` + `forkedAtHash`) plus the user's
   * chosen fork title. View choices are preserved — the fork inherits the
   * user's current view. No-op when no draft is active.
   */
  setDraftForkOf: (p: { newTitle: string; sourceId: number; sourceHash: string }) => void;
  addMember: (exposureId: number, qc: QueryClient) => void;
  removeMember: (index: number) => void;
  updateMember: (index: number, partial: Partial<DraftMember>) => void;
  reorderMembers: (newOrder: number[]) => void;
  resizeBands: (memberIdx: number, deltaPx: number, totalHeightPx: number) => void;
  resetBandHeights: () => void;
  /**
   * Cycle one peak's display state on a draft member (Plan §Phase 8.1):
   *   shown → labeled → hidden → shown (regular click)
   *   any   → hidden                  (alt+click)
   *
   * No-op when there's no active draft or `memberIdx` is out of range.
   */
  cyclePeakDisplayForMember: (memberIdx: number, peakId: number, altKey: boolean) => void;
  discardDraft: () => void;

  // Compare-page Phase 9 review-mode UI actions
  /**
   * Set the grouping mode on the active draft (C-4). Creates an empty draft
   * if none is active so the viewer can toggle without entering full edit mode
   * (spec §6.4 viewer escape hatch). effectiveGroupingMode(draft, comparison)
   * then surfaces the value to consumers.
   */
  setDraftViewGroupingMode: (mode: ActiveDraft["viewGroupingMode"]) => void;
  setShowPeakTicks: (show: boolean) => void;
  setShowPeakLabels: (show: boolean) => void;
  setHighlightedCompareMemberId: (id: number | undefined) => void;

  // Phase 12 — conflict modal slot
  setPendingConflict: (conflict: ConflictError | null) => void;

  // Permalink URL handling actions (spec §4.4 + §6)
  setStaleUrlContext: (ctx: StaleUrlContext | null) => void;
  setResolving: (v: boolean) => void;
  recoverFromStaleUrl: (opts: RecoverOpts) => void;
  /**
   * Atomic commit of a `/api/resolve` 200 response. Single setState so
   * `useUrlFromState` recomputes once — no cascading partial URL emits.
   * Arms `emitReplaceNext()` so the resulting state→URL emit is replace.
   */
  setResolveSuccess: (slots: ResolveSuccessSlots) => void;
  /** Mark the URL as an unknown frontend path (renders StaleUrlPage). */
  setStaleUnknownPath: (raw: string) => void;
  /** Atomic commit of a `/api/resolve` 404 response. Renders StaleUrlPage. */
  setStaleNotFound: (
    ctx: Extract<StaleUrlContext, { kind: "not_found" }>,
  ) => void;
}

/**
 * Wrap a state mutator so that any change to `activeDraft` is mirrored to
 * sessionStorage. We don't use Zustand's `persist` middleware for the draft
 * because we want sessionStorage (not localStorage) AND a separate schema
 * version, both of which `persist` can't accommodate alongside the LS_KEY
 * partition.
 */
function withDraftMirror(
  set: (partial: Partial<AppState>) => void,
  get: () => AppState,
) {
  return (next: ActiveDraftSlot): void => {
    set({ activeDraft: next });
    persistDraftToSession(next);
    void get; // unused, kept for symmetry with potential future read-paths
  };
}

export const useAppState = create<AppState>()(
  persist(
    (set, get) => {
      const setDraft = withDraftMirror(set, get);
      return {
        username: undefined,
        firstName: undefined,
        lastName: undefined,
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
        activePage: "index",
        tutorialSeen: false,
        theme: "dark",

        hoveredIndexId: undefined,
        hoveredPeakId: undefined,
        navModalOpen: false,
        navModalStep: "experiment",
        speculativeBuilder: null,
        compareXDomains: {},
        // Rehydrate the draft from sessionStorage at module-init time so
        // a tab reload restores edit-in-progress.
        activeDraft: loadDraftFromSession(),

        // Phase 9 — review-mode UI defaults. All per-tab; not persisted.
        showPeakTicks: true,
        showPeakLabels: true,
        highlightedCompareMemberId: undefined,

        // Phase 12 — conflict modal closed by default.
        pendingConflict: null,

        // Permalink URL handling — both ephemeral, default empty.
        staleUrlContext: null,
        resolving: false,

        setUsername: (username) => set({ username }),
        setUser: ({ username, firstName, lastName }) =>
          set({ username, firstName, lastName }),
        setActiveExperiment: (activeExperimentId) =>
          set({
            activeExperimentId,
            activeSampleId: undefined,
            activeExposureId: undefined,
            staleUrlContext: null,
          }),
        setActiveSample: (activeSampleId) =>
          set({ activeSampleId, activeExposureId: undefined, staleUrlContext: null }),
        setActiveExposure: (activeExposureId) => {
          // Only arm replace when (a) the value actually changes and (b) we're
          // on the Inspect page where the URL includes `?exposure=`. On Index,
          // the exposure isn't in the URL — PlotCard's auto-pick effect would
          // otherwise clobber the PUSH mode armed by an in-flight sample switch
          // (`,`/`.` keypress), breaking back-button history (issue #118).
          const s = get();
          if (s.activeExposureId !== activeExposureId && s.activePage === "inspect") {
            emitReplaceNext();
          }
          set({ activeExposureId, staleUrlContext: null });
        },
        setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
        setHoveredPeak: (hoveredPeakId) => set({ hoveredPeakId }),
        setActivePage: (activePage) => set({ activePage, staleUrlContext: null }),
        setTutorialSeen: (tutorialSeen) => set({ tutorialSeen }),
        setTheme: (theme) => set({ theme }),
        openNavModal: (step) =>
          set(step ? { navModalOpen: true, navModalStep: step } : { navModalOpen: true }),
        closeNavModal: () => set({ navModalOpen: false }),
        setNavModalStep: (navModalStep) => set({ navModalStep }),
        clearUsername: () => set({ username: undefined, firstName: undefined, lastName: undefined }),
        openSpeculativeBuilder: (exposureId) =>
          set({ speculativeBuilder: { exposureId } }),
        closeSpeculativeBuilder: () => set({ speculativeBuilder: null }),
        setCompareXDomain: (id, d) =>
          set({ compareXDomains: { ...get().compareXDomains, [id]: d } }),

        // ── Compare-draft actions ──────────────────────────────────────
        // Guard: re-calling on an already-new draft (id undefined) is a
        // no-op so the ComparePageEdit hydration effect can re-run without
        // clobbering an in-progress draft. Keeps the effect's deps array
        // exhaustive (no `draft` read inside the effect → no eslint-disable).
        startNewDraft: () => {
          const cur = get().activeDraft;
          if (cur !== null && cur.id === undefined) return;
          setDraft(emptyDraft());
        },
        startForkDraft: (comparison, qc) =>
          setDraft(fromComparisonAsFork(comparison, qc)),
        // Guard: re-loading the same comparison id is a no-op (don't clobber
        // an in-progress edit) — see startNewDraft above for rationale.
        loadDraftFromComparison: (comparison, qc) => {
          const cur = get().activeDraft;
          if (cur !== null && cur.id === comparison.id) return;
          setDraft(fromComparison(comparison, qc));
        },
        setDraftTitle: (title) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          setDraft({ ...cur, title });
        },
        setDraftDescription: (description) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          setDraft({ ...cur, description });
        },
        setDraftForkOf: (p) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          // Clearing `id` + `baseHash` flips the next submit onto the create
          // path (POST /api/comparisons, no expected_content_hash). The
          // `...cur` spread preserves members and view choices so the fork
          // inherits the user's current state.
          setDraft({
            ...cur,
            id: undefined,
            baseHash: undefined,
            forkedFromId: p.sourceId,
            forkedAtHash: p.sourceHash,
            title: p.newTitle,
          });
        },
        addMember: (exposureId, qc) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          const next: ActiveDraft = {
            ...cur,
            members: [
              ...cur.members,
              memberFromNewExposure(exposureId, cur.members.length, qc),
            ],
          };
          setDraft(next);
        },
        removeMember: (index) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          const filtered = cur.members.filter((_, i) => i !== index);
          // Renumber display_order so the contiguous range stays intact.
          const renumbered = filtered.map((m, i) => ({ ...m, display_order: i }));
          setDraft({ ...cur, members: renumbered });
        },
        updateMember: (index, partial) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          const members = cur.members.slice();
          if (index < 0 || index >= members.length) return;
          members[index] = { ...members[index]!, ...partial };
          setDraft({ ...cur, members });
        },
        reorderMembers: (newOrder) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          if (newOrder.length !== cur.members.length) return;
          const seen = new Set<number>();
          for (const idx of newOrder) {
            if (idx < 0 || idx >= cur.members.length || seen.has(idx)) return;
            seen.add(idx);
          }
          const reordered = newOrder.map((idx, i) => ({
            ...cur.members[idx]!,
            display_order: i,
          }));
          setDraft({ ...cur, members: reordered });
        },
        resizeBands: (memberIdx, deltaPx, totalHeightPx) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          if (memberIdx < 0 || memberIdx >= cur.members.length - 1) return;
          if (totalHeightPx <= 0) return;
          // Convert pixel delta to a band-height ratio. The dragged member
          // grows by deltaRatio = deltaPx/totalHeightPx; the next neighbor
          // shrinks by the same amount, preserving total band height.
          // Floors at a small minimum so a band can't collapse to zero.
          const deltaRatio = deltaPx / totalHeightPx;
          const MIN_HEIGHT = 0.1;
          const a = cur.members[memberIdx]!;
          const b = cur.members[memberIdx + 1]!;
          const newA = Math.max(MIN_HEIGHT, a.band_height + deltaRatio);
          const newB = Math.max(MIN_HEIGHT, b.band_height - deltaRatio);
          // If clamping ate part of the delta, propagate the actual change so
          // total height (sum of band_heights) stays approximately stable.
          const actualDelta = newA - a.band_height;
          const adjustedB = Math.max(MIN_HEIGHT, b.band_height - actualDelta);
          const members = cur.members.slice();
          members[memberIdx] = { ...a, band_height: newA };
          members[memberIdx + 1] = { ...b, band_height: adjustedB === newB ? newB : adjustedB };
          setDraft({ ...cur, members });
        },
        resetBandHeights: () => {
          const cur = get().activeDraft;
          if (cur === null) return;
          const members = cur.members.map((m) => ({ ...m, band_height: 1 }));
          setDraft({ ...cur, members });
        },
        cyclePeakDisplayForMember: (memberIdx, peakId, altKey) => {
          const cur = get().activeDraft;
          if (cur === null) return;
          if (memberIdx < 0 || memberIdx >= cur.members.length) return;
          const target = cur.members[memberIdx]!;
          const next = cyclePeakDisplay(target.peak_display, peakId, altKey);
          const members = cur.members.slice();
          members[memberIdx] = { ...target, peak_display: next };
          setDraft({ ...cur, members });
        },
        discardDraft: () => setDraft(null),

        // Phase 9 / C-4 — view-choice actions
        setDraftViewGroupingMode: (mode) => {
          const cur = get().activeDraft;
          // Viewer escape hatch (spec §6.4): if no draft is active, create an
          // empty one so the grouping preference can be carried without forcing
          // the user into full edit mode.
          const base = cur ?? emptyDraft();
          setDraft({ ...base, viewGroupingMode: mode });
        },
        setShowPeakTicks: (showPeakTicks) => set({ showPeakTicks }),
        setShowPeakLabels: (showPeakLabels) => set({ showPeakLabels }),
        setHighlightedCompareMemberId: (highlightedCompareMemberId) =>
          set({ highlightedCompareMemberId }),

        // Phase 12 — replace, never stack. A second 409 mid-modal updates
        // `current_state` in-place; the modal stays open.
        setPendingConflict: (pendingConflict) => set({ pendingConflict }),

        // Permalink URL handling actions (spec §4.4 + §6).
        // `recoverFromStaleUrl` is atomic: clears stale + sets active ids +
        // opens nav modal in one render-cycle commit so consumers don't see
        // an intermediate state.
        setStaleUrlContext: (staleUrlContext) => set({ staleUrlContext }),
        setResolving: (resolving) => set({ resolving }),
        recoverFromStaleUrl: (opts) => {
          emitReplaceNext();
          set((s) => ({
            staleUrlContext: null,
            activeExperimentId: opts.experimentId ?? s.activeExperimentId,
            activeSampleId: opts.sampleId ?? undefined,
            activeExposureId: undefined,
            navModalOpen: opts.openModal ?? true,
            navModalStep: opts.step,
          }));
        },
        setResolveSuccess: ({ page, experimentId, sampleId, exposureId }) => {
          emitReplaceNext();
          set({
            activePage: page,
            activeExperimentId: experimentId,
            activeSampleId: sampleId,
            activeExposureId: exposureId,
            staleUrlContext: null,
            resolving: false,
          });
        },
        setStaleUnknownPath: (raw) =>
          set({
            staleUrlContext: { kind: "unknown_path", raw },
            resolving: false,
          }),
        setStaleNotFound: (ctx) =>
          set({
            staleUrlContext: ctx,
            resolving: false,
          }),
      };
    },
    {
      name: LS_KEY,
      version: 3,
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        activePage: s.activePage,
        tutorialSeen: s.tutorialSeen,
        theme: s.theme,
      }),
    },
  ),
);
