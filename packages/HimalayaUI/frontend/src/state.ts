import { create } from "zustand";
import { persist } from "zustand/middleware";
import type { QueryClient } from "@tanstack/react-query";
import type { Comparison } from "./api";
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
import type { GroupingMode } from "./lib/comparison/coloring";

export const LS_KEY = "himalaya-ui:state";

export type PageId = "index" | "compare" | "inspect";
export type ThemeId = "dark" | "light";
export type NavModalStep = "experiment" | "sample";

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
   * Compare-page review/edit-mode grouping mode for trace coloring (Plan
   * §Phase 9, Task 9.2). Per-tab viewing preference — NOT persisted on the
   * comparison and NOT mirrored to sessionStorage. Default `"bySample"`.
   */
  groupingMode: GroupingMode;

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
  setGroupingMode: (mode: GroupingMode) => void;
  setShowPeakTicks: (show: boolean) => void;
  setShowPeakLabels: (show: boolean) => void;
  setHighlightedCompareMemberId: (id: number | undefined) => void;
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
        groupingMode: "bySample",
        showPeakTicks: true,
        showPeakLabels: true,
        highlightedCompareMemberId: undefined,

        setUsername: (username) => set({ username }),
        setUser: ({ username, firstName, lastName }) =>
          set({ username, firstName, lastName }),
        setActiveExperiment: (activeExperimentId) =>
          set({
            activeExperimentId,
            activeSampleId: undefined,
            activeExposureId: undefined,
          }),
        setActiveSample: (activeSampleId) =>
          set({ activeSampleId, activeExposureId: undefined }),
        setActiveExposure: (activeExposureId) => set({ activeExposureId }),
        setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
        setHoveredPeak: (hoveredPeakId) => set({ hoveredPeakId }),
        setActivePage: (activePage) => set({ activePage }),
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
        startNewDraft: () => setDraft(emptyDraft()),
        startForkDraft: (comparison, qc) =>
          setDraft(fromComparisonAsFork(comparison, qc)),
        loadDraftFromComparison: (comparison, qc) =>
          setDraft(fromComparison(comparison, qc)),
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

        // Phase 9 — review-mode UI actions
        setGroupingMode: (groupingMode) => set({ groupingMode }),
        setShowPeakTicks: (showPeakTicks) => set({ showPeakTicks }),
        setShowPeakLabels: (showPeakLabels) => set({ showPeakLabels }),
        setHighlightedCompareMemberId: (highlightedCompareMemberId) =>
          set({ highlightedCompareMemberId }),
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
