import { create } from "zustand";
import { persist } from "zustand/middleware";
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
import type { OrderRule, Series, IngestStatus } from "./api";

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

export type NavModalStep = "experiment" | "sample";

export type IngestProgressStatus = Exclude<IngestStatus, "idle">;
export interface IngestProgress {
  processed: number;
  total: number;
  status: IngestProgressStatus;
  /** Which segment of the progress bar is reporting ("discovery" / "analyzing" /
   *  "thumbnails"). SEPARATE from `status`, which selects the whole surface --
   *  see lib/ingestStages.ts. Undefined on an `ingest_started` frame, or from a
   *  backend that predates stage reporting; the UI then falls back to the plain
   *  single-track bar. */
  // `| undefined` is required, not redundant: exactOptionalPropertyTypes is on,
  // so `stage?: string` would reject an explicit `stage: undefined` — which is
  // exactly what the SSE listener passes when a frame carries no stage.
  stage?: string | undefined;
}

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
   *  state (mirrors the `staleUrlContext` non-persist rationale below). */
  hoveredQ: number | undefined;
  /** Plan D-7: ephemeral hypothetical-candidate preview. When set, the plot
   *  (combs ghost row + trace highlight) previews what the cart WOULD become
   *  if this candidate index were added — pure plot-only, NEVER a mutator,
   *  never an SSE event. Omitted from partialize (a momentary tab-local cue);
   *  cleared on hover-leave/blur so a stale ghost never masks the real cart. */
  previewIndexId: number | undefined;
  /** Per-experiment live-ingest progress (spec §9.3/§9.6). Ephemeral — written
   *  by the App.tsx SSE listener on `ingest_*` frames, read by the experiment
   *  header + LiveIngestUnfold. NOT persisted (omitted from partialize); a
   *  reload re-derives terminal state from the experiment's `ingest_status`. */
  ingestInFlight: Record<number, IngestProgress> | null;
  navModalOpen: boolean;
  navModalStep: NavModalStep;
  helpOverlayOpen: boolean;

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
   * Permalink URL handling slot (spec §4.4 + §6).
   * Ephemeral — not persisted. Non-null when the current URL points to a slug
   * that doesn't resolve (404 from `/api/resolve` or unknown path).
   */
  staleUrlContext: StaleUrlContext | null;

  // setters
  setUsername: (name: string) => void;
  setUser: (u: { username: string; firstName?: string | undefined; lastName?: string | undefined }) => void;
  setActiveExperiment: (id: number | undefined) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
  setHoveredIndex: (id: number | undefined) => void;
  setHoveredPeak: (id: number | undefined) => void;
  setHoveredQ: (q: number | undefined) => void;
  setPreviewIndex: (id: number | undefined) => void;
  setIngestProgress: (experimentId: number, progress: IngestProgress) => void;
  clearIngestProgress: (experimentId: number) => void;
  setTutorialSeen: (seen: boolean) => void;
  openNavModal: (step?: NavModalStep) => void;
  closeNavModal: () => void;
  openHelpOverlay: () => void;
  closeHelpOverlay: () => void;
  setNavModalStep: (step: NavModalStep) => void;
  clearUsername: () => void;

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
  /**
   * Replace the whole draft slot (or clear it with `null`). The primitive the
   * builder's undo/redo restores a prior snapshot through — restoring `null`
   * steps the first edit back out to read mode.
   */
  restoreSeriesDraft: (slot: SeriesDraftSlot) => void;

  // Compare-page Phase 9 review-mode UI actions (annotation toggles only)
  setShowPeakTicks: (show: boolean) => void;
  setShowPeakLabels: (show: boolean) => void;

  // Permalink URL handling actions (spec §4.4 + §6)
  setStaleUrlContext: (ctx: StaleUrlContext | null) => void;
  recoverFromStaleUrl: (opts: RecoverOpts) => void;
  /** Mark the URL as an unknown frontend path (renders StaleUrlPage). */
  setStaleUnknownPath: (raw: string) => void;
  /** Atomic commit of a `/api/resolve` 404 response. Renders StaleUrlPage. */
  setStaleNotFound: (
    ctx: Extract<StaleUrlContext, { kind: "not_found" }>,
  ) => void;
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
      const setSeriesDraft = withSeriesDraftMirror(set);
      return {
        username: undefined,
        firstName: undefined,
        lastName: undefined,
        activeExperimentId: undefined,
        activeSampleId: undefined,
        activeExposureId: undefined,
        tutorialSeen: false,

        hoveredIndexId: undefined,
        hoveredPeakId: undefined,
        hoveredQ: undefined,
        previewIndexId: undefined,
        ingestInFlight: null,
        navModalOpen: false,
        navModalStep: "experiment",
        helpOverlayOpen: false,
        // I3.5b — rehydrate the series-builder draft from sessionStorage.
        seriesDraft: loadSeriesDraftFromSession(),

        // Phase 9 — review-mode UI defaults. All per-tab; not persisted.
        showPeakTicks: true,
        showPeakLabels: true,

        // Permalink URL handling — ephemeral, default empty.
        staleUrlContext: null,

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
          set({
            activeSampleId,
            activeExposureId: undefined,
            staleUrlContext: null,
          }),
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
        setPreviewIndex: (previewIndexId) => set({ previewIndexId }),
        setIngestProgress: (experimentId, progress) =>
          set((s) => ({
            ingestInFlight: { ...(s.ingestInFlight ?? {}), [experimentId]: progress },
          })),
        clearIngestProgress: (experimentId) =>
          set((s) => {
            if (s.ingestInFlight === null) return {};
            const next = { ...s.ingestInFlight };
            delete next[experimentId];
            return { ingestInFlight: Object.keys(next).length === 0 ? null : next };
          }),
        setTutorialSeen: (tutorialSeen) => set({ tutorialSeen }),
        openNavModal: (step) =>
          set(step ? { navModalOpen: true, navModalStep: step } : { navModalOpen: true }),
        closeNavModal: () => set({ navModalOpen: false }),
        openHelpOverlay: () => set({ helpOverlayOpen: true }),
        closeHelpOverlay: () => set({ helpOverlayOpen: false }),
        setNavModalStep: (navModalStep) => set({ navModalStep }),
        clearUsername: () => set({ username: undefined, firstName: undefined, lastName: undefined }),

        // ── Series-builder draft actions (I3.5b) ─────────────────────────
        startSeriesDraftFromSeries: (series) => {
          const cur = get().seriesDraft;
          // Idempotent: keep an in-progress edit for the same series id.
          if (cur !== null && cur.id === series.id) return;
          setSeriesDraft(fromSeries(series));
        },
        discardSeriesDraft: () => setSeriesDraft(null),
        restoreSeriesDraft: (slot) => setSeriesDraft(slot),
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

        // Phase 9 — annotation toggle actions
        setShowPeakTicks: (showPeakTicks) => set({ showPeakTicks }),
        setShowPeakLabels: (showPeakLabels) => set({ showPeakLabels }),

        // Permalink URL handling actions (spec §4.4 + §6).
        // `recoverFromStaleUrl` is atomic: clears stale + sets active ids +
        // opens nav modal in one render-cycle commit so consumers don't see
        // an intermediate state.
        setStaleUrlContext: (staleUrlContext) => set({ staleUrlContext }),
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
          set({ staleUrlContext: { kind: "unknown_path", raw } }),
        setStaleNotFound: (ctx) =>
          set({ staleUrlContext: ctx }),
      };
    },
    {
      name: LS_KEY,
      // I5.1 (#182) dropped `activePage` from partialize (its field is gone).
      // I5.2 (#183) bumped `version` 3 → 4 WITH the `migrate` below, which
      // formally strips any lingering `activePage` key from a pre-cutover blob.
      // R0a (#221): "The Print" is the single identity — the dark `theme` field
      // is gone, so it leaves partialize and `version` bumps 4 → 5; the migrate
      // now ALSO strips a lingering `theme` key. A version bump WITHOUT a
      // migrate would make Zustand discard the whole persisted blob — the wipe
      // risk; the migrate preserves every surviving pref
      // (username/tutorialSeen/active*/…).
      version: 5,
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        tutorialSeen: s.tutorialSeen,
      }),
      // I5.2 (#183): runs BEFORE `merge` on rehydrate (zustand v4 order is
      // migrate → merge). Returns persisted DATA only — `merge`
      // (mergePersistedState) re-attaches the store actions afterward via
      // `...current`. The transforms across every prior version are "drop the
      // dead `activePage` key" (I5.2 #183) and "drop the retired `theme` key"
      // (R0a #221) — both are unconditional key-strips, so there is no
      // `switch (version)`.
      migrate: (persisted, _version) => {
        // WIPE-GUARD: a non-object / malformed blob is returned UNTOUCHED —
        // never `{}`/`undefined`. Handing `{}` here would be the only way this
        // migrate could itself partial-wipe prefs; returning the original lets
        // `merge` fold whatever survived (or fall back to defaults).
        if (persisted && typeof persisted === "object") {
          const {
            activePage: _activePage,
            theme: _theme,
            ...rest
          } = persisted as Record<string, unknown>;
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
