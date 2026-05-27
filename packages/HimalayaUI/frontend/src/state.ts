import { create } from "zustand";
import { persist } from "zustand/middleware";
import type { ConflictError } from "./api";
import {
  emptyDraft,
  loadDraftFromSession,
  persistDraftToSession,
  type ActiveDraft,
  type ActiveDraftSlot,
  type DraftMember,
} from "./lib/comparison/draft";
import {
  loadSeriesDraftFromSession,
  persistSeriesDraftToSession,
  type SeriesDraftSlot,
} from "./lib/series/seriesDraft";
import {
  fromSeries,
  addSampleToRecipe,
  removeRecipeRow,
  reorderRecipe,
} from "./lib/series/seriesDraftFactories";
import type { OrderRule, Series } from "./api";

export const LS_KEY = "himalaya-ui:state";

// I5.1 (#182): the dual-nav `activePage` model — `PageId`, `VALID_PAGE_IDS`,
// `coerceActivePage`, the `activePage` store field + `setActivePage` — is
// deleted. (The resolve setter's former `page` slot went with it; I5.3 (#184)
// then removed the now-unreachable `setResolveSuccess` action entirely.) Every legacy AppShell
// surface (Inspect #163, Index #181, Compare #177) is retired and the app is a
// single URL-routed shell, so nothing renders off `activePage` anymore. The
// `activePage: s.activePage` line is also removed from `partialize` (forced by
// the field deletion). I5.2 (#183) then bumped the persist `version` 3 → 4 with
// the `migrate` below, which formally strips any lingering `activePage` key from
// a pre-cutover blob while preserving every surviving pref.

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

export interface AppState {
  // persisted
  username: string | undefined;
  firstName: string | undefined;
  lastName: string | undefined;
  activeExperimentId: number | undefined;
  activeSampleId: number | undefined;
  activeExposureId: number | undefined;
  tutorialSeen: boolean;
  theme: ThemeId;

  // ephemeral (not persisted to localStorage)
  hoveredIndexId: number | undefined;
  hoveredPeakId: number | undefined;
  /** Ephemeral q-link hover (focus workspace, #180). The hovered q-value that
   *  cross-highlights the trace peak and the detector ring.
   *
   *  NOT persisted — deliberately omitted from `partialize`, like the other
   *  hovered* fields; a hovered q is a momentary, tab-local UI cue that would
   *  be meaningless to replay across reloads. It is also invisible to the
   *  mutation queue and never reaches a server payload — pure client hover
   *  state (mirrors the `pendingConflict` non-persist rationale below). */
  hoveredQ: number | undefined;
  navModalOpen: boolean;
  navModalStep: NavModalStep;
  // Speculative builder: when non-null, the modal is open for this exposure.
  // All builder form state (phase, anchor peak, ratio) is local to the
  // SpeculativeBuilder component — only the open/close gate lives in store
  // because PhasePanel needs to mount/unmount the modal.
  speculativeBuilder: { exposureId: number } | null;

  // ── Compare-era draft / view slice — KEPT (I3.6 #177 deviation; see below) ──
  //
  // The I3.6 plan (§3.2 / §8) resolved to REMOVE this slice when the Compare
  // page was retired. That resolution was written against the PRE-I3.5b tree;
  // it is unsafe against the tree this PR actually rebases onto. I3.5b built the
  // series builder ON TOP OF this slice, so a chunk of it is now live, shared
  // infrastructure — NOT dead Compare-only state. grep-verified surviving
  // (non-test, non-Compare) consumers:
  //   - `showPeakTicks` / `showPeakLabels` (+ setters): read directly by
  //     `SeriesBuilderPage.tsx`, `AnnotationToggles`, `MultiTracePlot`,
  //     `MemberTraceLayer`, and the figure-export adapters/marks.
  //   - `compareXDomains`: read by `SeriesBuilderPage.tsx`.
  //   - `activeDraft` + `updateMember` / `reorderMembers` / `resizeBands` /
  //     `setDraftViewGroupingMode` / `highlightedCompareMemberId` (+ setter):
  //     read by the shared render components the series builder mounts
  //     (`MemberMetaRow`, `MemberMetaGutter`, `BandResizeDivider`,
  //     `GroupingModeToggle`).
  // A genuinely-dead SUBSET remains (the create/fork/membership actions only the
  // deleted Compare page drove: `startNewDraft`, `startForkDraft`,
  // `loadDraftFromComparison`, `setDraftForkOf`, `addMember`, `removeMember`,
  // `discardDraft`, `setCompareXDomain`, `resetBandHeights`,
  // `cyclePeakDisplayForMember`, …). Pruning only that subset is interconnected
  // and type-shared (`ActiveDraft`) and risks the series builder's draft-backed
  // editing; it is DELIBERATELY DEFERRED to I5.3's (#184) dead-code sweep, which
  // runs against a stable post-cutover tree. I5.1 (#182) deleted the `activePage`
  // model but leaves this draft slice intact. (See the PR's coordination note.)
  /**
   * Compare-page draft slot (Plan §Phase 4, Task 4.3). Single slot — only
   * one comparison can be in edit mode at a time per tab. Mirrored to
   * sessionStorage with a schema version (see `lib/comparison/draft.ts`).
   * Tab close loses the draft, which is acceptable for v1 per the spec.
   */
  activeDraft: ActiveDraftSlot;

  /**
   * Series-builder draft slot (I3.5b). A SEPARATE namespace from `activeDraft`
   * — a series recipe edits `series_samples` membership, a different shape from
   * a comparison's plate, and overloading one slot would let the two flows
   * clobber each other across tabs. Single slot (one series in edit mode per
   * tab); mirrored to sessionStorage with its own schema key (see
   * `lib/series/seriesDraft.ts`); NOT persisted in localStorage / partialize.
   */
  seriesDraft: SeriesDraftSlot;

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
  setHoveredQ: (q: number | undefined) => void;
  setTutorialSeen: (seen: boolean) => void;
  setTheme: (theme: ThemeId) => void;
  openNavModal: (step?: NavModalStep) => void;
  closeNavModal: () => void;
  setNavModalStep: (step: NavModalStep) => void;
  clearUsername: () => void;
  openSpeculativeBuilder: (exposureId: number) => void;
  closeSpeculativeBuilder: () => void;

  // Compare-era draft actions — KEPT (consumed by the shared series-builder
  // render core: MemberMetaRow / MemberMetaGutter / BandResizeDivider). The
  // Compare-only create/fork/membership sub-actions were removed in I5.3 (#184).
  updateMember: (index: number, partial: Partial<DraftMember>) => void;
  reorderMembers: (newOrder: number[]) => void;
  resizeBands: (memberIdx: number, deltaPx: number, totalHeightPx: number) => void;

  // ── Series-builder draft actions (I3.5b) ───────────────────────────────
  /**
   * Seed the series draft from a loaded series. Idempotent: a no-op when a
   * draft for the same series id is already active (so a hydration effect can
   * re-run without clobbering an in-progress edit).
   */
  startSeriesDraftFromSeries: (series: Series) => void;
  discardSeriesDraft: () => void;
  setSeriesDraftTitle: (title: string) => void;
  setSeriesDraftDescription: (description: string) => void;
  setSeriesOrderingVariable: (value: string | null) => void;
  setSeriesOrderRule: (rule: OrderRule) => void;
  /** Append a sample to the recipe (negative placeholder id). No-op if no draft. */
  addSeriesSample: (sampleId: number) => void;
  /** Remove a recipe row by its local id. No-op if no draft. */
  removeSeriesSample: (rowId: number) => void;
  /** Move a recipe row from index `from` to index `to`. No-op if no draft. */
  reorderSeriesSample: (from: number, to: number) => void;

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

/**
 * Series-draft equivalent of `withDraftMirror` (I3.5b) — mirrors every
 * `seriesDraft` change to sessionStorage under its own schema key, for the
 * same reasons (sessionStorage + separate version).
 */
function withSeriesDraftMirror(set: (partial: Partial<AppState>) => void) {
  return (next: SeriesDraftSlot): void => {
    set({ seriesDraft: next });
    persistSeriesDraftToSession(next);
  };
}

export const useAppState = create<AppState>()(
  persist(
    (set, get) => {
      const setDraft = withDraftMirror(set, get);
      const setSeriesDraft = withSeriesDraftMirror(set);
      return {
        username: undefined,
        firstName: undefined,
        lastName: undefined,
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
        tutorialSeen: false,
        theme: "dark",

        hoveredIndexId: undefined,
        hoveredPeakId: undefined,
        hoveredQ: undefined,
        navModalOpen: false,
        navModalStep: "experiment",
        speculativeBuilder: null,
        // Rehydrate the draft from sessionStorage at module-init time so
        // a tab reload restores edit-in-progress.
        activeDraft: loadDraftFromSession(),
        // I3.5b — same rehydration for the series-builder draft.
        seriesDraft: loadSeriesDraftFromSession(),

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
          // Inspect — the only surface that put the exposure in the URL via
          // `?exposure=` — is retired (#163). The replace-mode arming guarded
          // on `activePage === "inspect"` (issue #118) goes with it; no live
          // surface carries the exposure in the URL, so just set the value.
          set({ activeExposureId, staleUrlContext: null });
        },
        setHoveredIndex: (hoveredIndexId) => set({ hoveredIndexId }),
        setHoveredPeak: (hoveredPeakId) => set({ hoveredPeakId }),
        setHoveredQ: (hoveredQ) => set({ hoveredQ }),
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

        // ── Compare-era draft actions — KEPT ───────────────────────────
        // Consumed by the shared series-builder render core (MemberMetaRow,
        // MemberMetaGutter, BandResizeDivider). The Compare-only create/fork/
        // membership sub-actions were removed in I5.3 (#184).
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

        // ── Series-builder draft actions (I3.5b) ─────────────────────────
        startSeriesDraftFromSeries: (series) => {
          const cur = get().seriesDraft;
          // Idempotent: keep an in-progress edit for the same series id.
          if (cur !== null && cur.id === series.id) return;
          setSeriesDraft(fromSeries(series));
        },
        discardSeriesDraft: () => setSeriesDraft(null),
        setSeriesDraftTitle: (title) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft({ ...cur, title });
        },
        setSeriesDraftDescription: (description) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft({ ...cur, description });
        },
        setSeriesOrderingVariable: (value) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft({ ...cur, orderingVariable: value });
        },
        setSeriesOrderRule: (rule) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft({ ...cur, orderRule: rule });
        },
        addSeriesSample: (sampleId) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft(addSampleToRecipe(cur, sampleId));
        },
        removeSeriesSample: (rowId) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft(removeRecipeRow(cur, rowId));
        },
        reorderSeriesSample: (from, to) => {
          const cur = get().seriesDraft;
          if (cur === null) return;
          setSeriesDraft(reorderRecipe(cur, from, to));
        },

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
          set((s) => ({
            staleUrlContext: null,
            activeExperimentId: opts.experimentId ?? s.activeExperimentId,
            activeSampleId: opts.sampleId ?? undefined,
            activeExposureId: undefined,
            navModalOpen: opts.openModal ?? true,
            navModalStep: opts.step,
          }));
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
      // I5.1 (#182) dropped `activePage` from partialize (its field is gone).
      // I5.2 (#183) bumps `version` 3 → 4 WITH the `migrate` below, which
      // formally strips any lingering `activePage` key from a pre-cutover blob.
      // A version bump WITHOUT a migrate would make Zustand discard the whole
      // persisted blob — the wipe risk; the migrate preserves every surviving
      // pref (username/theme/tutorialSeen/…).
      version: 4,
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        tutorialSeen: s.tutorialSeen,
        theme: s.theme,
      }),
      // I5.2 (#183): runs BEFORE `merge` on rehydrate (zustand v4 order is
      // migrate → merge). Returns persisted DATA only — `merge`
      // (mergePersistedState) re-attaches the store actions afterward via
      // `...current`. The only transform across every prior version is "drop
      // the dead `activePage` key", so there is no `switch (version)`.
      migrate: (persisted, _version) => {
        // WIPE-GUARD: a non-object / malformed blob is returned UNTOUCHED —
        // never `{}`/`undefined`. Handing `{}` here would be the only way this
        // migrate could itself partial-wipe prefs; returning the original lets
        // `merge` fold whatever survived (or fall back to defaults).
        if (persisted && typeof persisted === "object") {
          const { activePage: _activePage, ...rest } = persisted as Record<
            string,
            unknown
          >;
          return rest as unknown as AppState;
        }
        return persisted as AppState;
      },
      merge: mergePersistedState,
    },
  ),
);

/** persist `merge` — replicates zustand's default shallow merge
 *  ({ ...current, ...persisted }). Kept (rather than relying on zustand's
 *  default) so the merge strategy is explicit and unit-testable. The former
 *  `activePage` coercion is gone with the dual-nav model (I5.1, #182). This
 *  merge runs AFTER `migrate` on rehydrate (zustand v4) and is unchanged by
 *  the I5.2 (#183) `version` 3 → 4 bump — the strip happens in `migrate`. */
export function mergePersistedState(
  persisted: unknown,
  current: AppState,
): AppState {
  return { ...current, ...(persisted as Partial<AppState> | undefined) };
}
