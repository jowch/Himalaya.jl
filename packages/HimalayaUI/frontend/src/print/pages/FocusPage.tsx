import { useMemo, useState, useCallback, useRef, useEffect } from "react";
import { useNavigate, useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { TracePlate } from "../components/TracePlate";
import { DetectorPanel } from "../components/DetectorPanel";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { CombsPanel, type CombView } from "../components/CombsPanel";
import type { PlotPeak, PeakFocusRequest } from "../plot/marks/PlotPeaks";
import { nextFocusPeakId } from "../plot/peakFocusOrder";
import { AssignmentRail } from "../components/AssignmentRail";
import { AssignmentCart } from "../components/AssignmentCart";
import { PhaseBlock } from "../components/PhaseBlock";
import { CandidateRow, CandidateList } from "../components/CandidateRow";
import { FormFactorRow } from "../components/FormFactorRow";
import { CustomIndexModal } from "../components/CustomIndexModal";
import { HintText, EmptyState, Button } from "../ui";
import { ExportButton } from "../components/ExportButton";
import { useFigureExport } from "../components/useFigureExport";
import { buildCleanFigureSvg, type FigureTraceKey } from "../export/cleanFigureSvg";
import { buildFocusFigureRow } from "../export/focusFigureRow";
import { phaseSegment } from "../export/seriesFigureKeys";
import { phaseHex } from "../export/traceColors";
import { waterfallQDomain } from "../waterfall/waterfallModel";
import {
  toTraceModel,
  peakClickAction,
  toDetectorRings,
  toCombSeries,
  buildDetectorCalibration,
  CUSTOM_SYMS,
  customIndexPreview,
} from "./focusAdapters";
import { phaseColor } from "../../phases";
import { buildExposureImageUrl, toGalleryExposures } from "./loupeAdapters";
import {
  useCorpusSamples,
  useExperiment,
  useExposures,
  useTrace,
  usePeaks,
  useIndices,
  useAssignment,
  useAddPeak,
  useRemovePeak,
  useSetPeakExcluded,
  useAddAssignmentPhase,
  useRemoveAssignmentPhase,
  useSetAssignmentState,
  useCommitCustomIndex,
} from "../../queries";
import { useAppState } from "../../state";
import { useSyncActiveSampleFromRoute } from "../../hooks/useSyncActiveSampleFromRoute";
import { useAutoPickExposure, noUsableExposureState, resolveActiveExposure } from "../../hooks/useAutoPickExposure";
import { useDocumentTitle } from "../../hooks/useDocumentTitle";
import { useExperimentSiblings } from "../../hooks/useExperimentSiblings";
import { useListCursor } from "../interaction/useListCursor";
import { useStepperOnly } from "../interaction/useStepperOnly";
import { usePageActions } from "../interaction/usePageActions";
import { core, page } from "../interaction/core";
import { deriveActiveIndices } from "../../lib/assignment";
import { sanitizeDashes } from "../../lib/copy";
import { basisFor } from "../../lib/customIndex";
import { seriesRatio, ratioTerm } from "../../lib/seriesRatio";
import { announce } from "../../lib/announce";
import { showToast } from "../../lib/toast";
import type { Trace, IndexEntry } from "../../api";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// Boneyard fixture — a real render with mock props so the headless capture CLI
// measures the greenfield Focus body. Uses static composites with no-op handlers
// and a small synthetic trace so the plot draws real geometry.
const FIXTURE_TRACE = {
  trace: {
    q: [0.02, 0.04, 0.06, 0.1, 0.15, 0.2],
    I: [12000, 8000, 4800, 2400, 1200, 600],
    sigma: [80, 55, 35, 20, 12, 7],
  },
  peaks: [] as PlotPeak[],
  phase: null,
};
// EX-SPACING: the focus page grids are skeleton↔real ALIGNMENT INVARIANTS — the
// loading placeholder (FOCUS_FIXTURE) must use the SAME templates as the loaded
// render or the load→loaded transition shifts. One source each (literal Tailwind
// class strings so the JIT scanner still sees them).
//   - work column + 350px Notes/phase rail
//   - the detector / combs two-up split below the trace plate
const FOCUS_PAGE_GRID = "grid grid-cols-1 lg:grid-cols-[minmax(0,1fr)_350px]";
const FOCUS_SPLIT_GRID =
  "grid grid-cols-1 gap-5 lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]";
const FOCUS_FIXTURE = (
  <div className={FOCUS_PAGE_GRID}>
    {/* work column */}
    <div className="min-w-0 px-8 pt-7 pb-13">
      <div className="mx-auto flex min-w-0 max-w-[1180px] flex-col gap-[18px]">
        <TracePlate
          kicker="Integration"
          title="Sample"
          subtitle="JC000 · Experiment · frame-001"
          trace={FIXTURE_TRACE}
          scale="log"
          onScaleChange={() => {}}
          onAutoFit={() => {}}
          onToggleAddPeak={() => {}}
        />
        <div className={FOCUS_SPLIT_GRID}>
          <DetectorPanel src={null} />
          <div className="hidden lg:flex min-h-0 flex-col">
            <CombsPanel
              assigned={[]}
              view="comb"
              onViewChange={() => {}}
            />
          </div>
        </div>
      </div>
    </div>

    {/* rail column = AssignmentRail (owns the rail shell) */}
    <AssignmentRail
      assignment={
        <AssignmentCart>
          <PhaseBlock phase="Pn3m" score={0.87} meta="a = 142 Å · 5 reflections" />
        </AssignmentCart>
      }
      candidates={
        <CandidateList>
          <CandidateRow phase="Pn3m" score={0.87} why="explains 5 peaks · in the call" selected />
          <CandidateRow phase="Im3m" score={0.61} why="explains 3 peaks" />
        </CandidateList>
      }
      candidatesNote="A sample can be multiphasic, so check every phase that fits."
    />
  </div>
);

/**
 * FocusPage (greenfield) — the Focus workspace at /sample/:sampleId.
 *
 * Mirrors the 2-column workspace grid (work · rail)
 * but assembled entirely from src/print composites + the carried data layer.
 * The route param seeds Zustand `activeSampleId` (useSyncActiveSampleFromRoute);
 * useAutoPickExposure seeds `activeExposureId`. All hooks are called
 * unconditionally (hooks rule) before any early return.
 */
export function FocusPage(): JSX.Element {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();

  // ── route → store seeding ──────────────────────────────────────────────────
  // The status guards against the mid-session lie: a bogus /sample/:id never
  // seeds the store, so activeSampleId would keep pointing at the previous
  // sample. routeStatus "unknown" forces the not-found branch instead.
  const routeStatus = useSyncActiveSampleFromRoute();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const storedExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  useAutoPickExposure(activeSampleId);

  // ── data ────────────────────────────────────────────────────────────────────
  const corpusQ = useCorpusSamples();
  const corpusSample =
    activeSampleId !== undefined
      ? corpusQ.data?.find((s) => s.id === activeSampleId)
      : undefined;
  const experimentId = corpusSample?.experiment_id;

  const experimentQ = useExperiment(experimentId ?? 0);
  const exposuresQ = useExposures(activeSampleId);

  // ── detector beam-center calibration ──────────────────────────────────────────
  // The raw detector pixel dims + display orientation are reported by
  // DetectorImage after the async canvas fetch resolves the image-route headers
  // (Task 5/6). Hold them here so the calibration recomputes when either the
  // experiment beamline params or the raw size change. `imageAspect` MUST derive
  // from the RAW header dims (not a decoded-bitmap size) so the overlay aspect is
  // the true detector aspect. Until a frame loads, rawSize is null → calibration
  // is null → the DetectorPanel centered fallback.
  const [rawSize, setRawSize] = useState<{ w: number; h: number } | null>(null);
  const [detOrient, setDetOrient] = useState<"portrait" | "landscape">("portrait");
  const handleRawSize = useCallback(
    (w: number, h: number) => setRawSize((p) => (p?.w === w && p?.h === h ? p : { w, h })),
    [],
  );
  const calibration = buildDetectorCalibration(experimentQ.data, rawSize);
  const imageAspect = rawSize ? rawSize.w / rawSize.h : undefined;

  // FO-NAV-SKELETON: resolve the active exposure at RENDER time, not just in
  // useAutoPickExposure's post-paint effect. A sample switch clears the stored
  // id (the setActiveSample cascade); reading it raw would leave a one-render
  // gap where activeExposureId is undefined even when the new sample's exposures
  // are already cached — and that gap trips the skeleton's transition fade on
  // cached back-and-forth. Resolving here means a cached navigation keys the
  // trace/peaks queries to their (warm) data in the same frame, so the skeleton
  // only ever shows when data is genuinely absent. The effect still runs to
  // persist this id to the store for the exposure stepper and deliberate
  // switches; render reads the derived value so it never lags the store.
  const activeExposureId = resolveActiveExposure(storedExposureId, exposuresQ.data);

  // Inter-sample order: siblings drives the URL-driven sampleStepper.
  const { siblings } = useExperimentSiblings();

  const traceQ = useTrace(activeExposureId);
  const peaksQ = usePeaks(activeExposureId);
  const indicesQ = useIndices(activeExposureId);
  const assignmentQ = useAssignment(activeExposureId);

  // ── mutators (all exposure-scoped) ───────────────────────────────────────────
  const addPeak = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);
  const setPeakExcluded = useSetPeakExcluded(activeExposureId ?? 0);
  const addAssignmentPhase = useAddAssignmentPhase(activeExposureId ?? 0);
  const removeAssignmentPhase = useRemoveAssignmentPhase(activeExposureId ?? 0);
  const setAssignmentState = useSetAssignmentState(activeExposureId ?? 0);
  const commitCustomIndex = useCommitCustomIndex(activeExposureId ?? 0);

  // ── page-owned state ─────────────────────────────────────────────────────────
  const [scale, setScale] = useState<"log" | "lin">("log");
  const [addArmed, setAddArmed] = useState(false);
  const [xDomain, setXDomain] = useState<[number, number] | null>(null);
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  const [combView, setCombView] = useState<CombView>("comb");
  // preview-vs-cursor split (resolutions §"preview-vs-cursor reconciliation"):
  // previewWasExplicit = user keyboard-navigated/clicked a candidate (sticky).
  // hoverPreviewId = transient mouse hover. previewIndexId is DERIVED below.
  const [previewWasExplicit, setPreviewWasExplicit] = useState(false);
  const [hoverPreviewId, setHoverPreviewId] = useState<number | undefined>(undefined);

  // FO-NAV-STATE: FocusPage is NOT remounted on a same-route [ / ] sample step
  // (React Router reuses the routed element for the new :sampleId), so
  // page-owned interaction state would survive a sample switch. Reset the
  // per-sample bits whenever the active sample changes: the "+ Peak" arm (else
  // the first click on the next sample's trace silently mutates ITS peaks), the
  // manual zoom window (else the next trace can render outside a stale x-domain),
  // and the candidate preview (else the stale previewIndexId — which renders
  // nothing, since the old index id never re-matches — still makes the Esc ladder
  // eat the first Escape as a no-op clear instead of backing out to the sheet).
  // scale / combView are sticky preferences and intentionally persist.
  useEffect(() => {
    setAddArmed(false);
    setXDomain(null);
    setPreviewWasExplicit(false);
    setHoverPreviewId(undefined);
  }, [activeSampleId]);

  // ── keyboard focus re-anchor after a destructive peak edit (WCAG 2.4.3) ──────
  // Removing a peak via the keyboard unmounts its mark; without this, focus
  // falls to <body>. On a remove we compute the surviving neighbour (q-order)
  // and bump a nonce-keyed request; TracePlate → TracePlot → PlotPeaks consumes
  // it in a layout effect AFTER the re-render and moves focus to that mark. The
  // nonce ref drives a NEW request each time so two removes onto the same
  // survivor still re-fire; a foreign/SSE re-render leaves the nonce alone, so
  // it never steals focus. Null survivor → focus the "+ Peak" button instead.
  const [focusRequest, setFocusRequest] = useState<PeakFocusRequest | undefined>(undefined);
  const focusNonce = useRef(0);
  const requestPeakFocus = useCallback((id: number | null) => {
    focusNonce.current += 1;
    setFocusRequest({ id, nonce: focusNonce.current });
  }, []);

  // custom-index modal
  const [customOpen, setCustomOpen] = useState(false);
  const [customSym, setCustomSym] = useState<string>(CUSTOM_SYMS[0]!.name);
  const [customParam, setCustomParam] = useState<string>(
    String(CUSTOM_SYMS[0]!.min),
  );

  // ── derived ──────────────────────────────────────────────────────────────────
  const peaks = peaksQ.data ?? [];
  const indices = indicesQ.data ?? [];
  const assignment = assignmentQ.data;
  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);
  const activeExposure = exposures.find((e) => e.id === activeExposureId);
  // The candidate cursor list (speculatives excluded — they are the custom-index
  // outputs, toggled not previewed). Defined here (above the not-found early
  // return) so hooks that follow can close over it.
  const candidatePool = indices.filter((i) => i.kind !== "speculative");

  const activeIndices = useMemo(
    () => deriveActiveIndices(assignment, indices),
    [assignment, indices],
  );
  const memberIds = useMemo(
    () => new Set(activeIndices.map((ix) => ix.id)),
    [activeIndices],
  );

  // ── candidate cursor (headless — no roving rows; focus lives on scope) ───────
  const candidatePoolIds = useMemo(
    () => candidatePool.map((i) => i.id),
    [candidatePool],
  );
  const toggleAssignmentForId = useCallback(
    (id: number) => {
      const ix = indices.find((i) => i.id === id);
      if (!ix) return;
      if (memberIds.has(ix.id)) {
        removeAssignmentPhase.mutate(ix.id);
        announce(`${ix.phase} removed from the call`);
      } else {
        addAssignmentPhase.mutate(ix.id);
        announce(`${ix.phase} added to the call`);
      }
    },
    [indices, memberIds, removeAssignmentPhase, addAssignmentPhase],
  );
  const candidateCursor = useListCursor({
    ids: candidatePoolIds,
    onActivate: toggleAssignmentForId,
    // Any candidate move — ←/→ OR the dock Candidate stepper — makes the comb
    // preview explicit (and enables Apply). Centralizing it here keeps the two
    // controls consistent (the stepper used to skip it).
    onMove: () => setPreviewWasExplicit(true),
    stepperLabel: "Candidate",
    stepperTestIdBase: "candidate",
    axis: "horizontal",
  });

  // Derived preview: hover wins (transient), else the sticky keyboard preview,
  // else nothing. Existing readers of previewIndexId keep working unchanged.
  const previewIndexId: number | undefined =
    hoverPreviewId ?? (previewWasExplicit ? (candidateCursor.cursorId ?? undefined) : undefined);

  // ── sample stepper (URL-driven, extra) ───────────────────────────────────────
  const sampleStepper = useStepperOnly({
    ids: siblings.map((s) => s.id),
    currentId: activeSampleId,
    onGo: (id) => navigate(`/sample/${id}`),
    label: "Sample",
    testIdBase: "sample",
    axis: "vertical",
  });

  // ── scope container (focus anchor) ──────────────────────────────────────────
  // scopeEl: imperative .focus() in escapeLadder (WCAG 2.4.3 re-anchor).
  // scopeRef: callback ref — ONLY stores the element; focus-on-attach is NOT
  //   used here because the scope div renders inside a <Skeleton> and attaches
  //   while visibility:hidden (boneyard overlay), causing Chromium to silently
  //   drop the .focus() call. Focus is instead deferred to a useEffect that
  //   fires once isLoading transitions to false (the skeleton reveals content).
  const scopeEl = useRef<HTMLDivElement | null>(null);
  const focusedSampleRef = useRef<number | null>(null);
  const scopeRef = useCallback((el: HTMLDivElement | null) => {
    scopeEl.current = el;
  }, []);

  // ── action declaration ───────────────────────────────────────────────────────
  const fromSeries = searchParams.get("from") === "series";
  const backLabel = fromSeries ? "Series" : "Corpus";
  const goBack = useCallback(() => {
    navigate(fromSeries ? "/series" : `/experiments/${experimentId}/corpus`);
  }, [navigate, fromSeries, experimentId]);

  // Escape ladder (innermost first): disarm addPeak → clear sticky preview → leave.
  const escapeLadder = useCallback(() => {
    if (addArmed) {
      setAddArmed(false);
      scopeEl.current?.focus({ preventScroll: true });
      return;
    }
    if (previewWasExplicit) {
      setPreviewWasExplicit(false);
      return;
    }
    goBack();
  }, [addArmed, previewWasExplicit, goBack]);

  usePageActions({
    cursor: candidateCursor,
    extraSteppers: [sampleStepper],
    // ←/→ move the candidate cursor; ↑/↓ step the sample. Scope-exempt so arrows
    // control the surface wherever focus sits, instead of scrolling the page.
    arrowHandler: (e) => {
      if (e.key === "ArrowLeft") { e.preventDefault(); candidateCursor.moveBy(-1); }   // moveBy → onMove sets previewWasExplicit
      else if (e.key === "ArrowRight") { e.preventDefault(); candidateCursor.moveBy(1); }
      else if (e.key === "ArrowUp") { e.preventDefault(); sampleStepper.onPrev(); }
      else if (e.key === "ArrowDown") { e.preventDefault(); sampleStepper.onNext(); }
    },
    actions: [
      core("back", { label: backLabel, run: escapeLadder, dock: true }),
      core("openFocus", {
        label: "Apply",
        run: () => candidateCursor.activate(),
        dock: "primary",
        enabled: () => candidateCursor.cursorId !== null && previewWasExplicit,
      }),
      core("openLoupe", {
        run: () => {
          if (activeSampleId !== undefined)
            navigate(`/sample/${activeSampleId}/loupe`);
        },
        dock: true,
        enabled: () => activeSampleId !== undefined,
      }),
      page("addPeak", {
        label: "+ Peak",
        keys: ["p"],
        group: "Edit",
        dock: true,
        run: () => setAddArmed((v) => !v),
      }),
    ],
  });

  // Primary assigned phase: the first active index when the assignment is
  // `indexed` (mirrors the legacy trace-colour derivation, which reads the
  // active set off the assignment cart). Null otherwise → ink-faint trace.
  const phase =
    assignment?.state === "indexed" && activeIndices.length > 0
      ? activeIndices[0]!.phase
      : null;

  // Candidate-hover preview: light up the peaks the hovered candidate CLAIMS.
  // Those peaks recolour to the candidate's phase hue (a preview of "pick this
  // and these become this phase") while every other peak dims, so you can read
  // at a glance which peaks a candidate's score is built on. Replaces the old
  // active-phase "losing peak" dim, which went dark once the durable assignment
  // stopped auto-seeding (no active phase → nothing to lose → nothing lit).
  const previewIx = indices.find((i) => i.id === previewIndexId);
  const previewClaim = useMemo(
    () => (previewIx ? new Set(previewIx.peaks.map((p) => p.peak_id)) : undefined),
    [previewIx],
  );

  const traceModel = useMemo(() => {
    const base = toTraceModel(traceQ.data ?? EMPTY_TRACE, peaks, phase);
    if (!previewIx || !previewClaim) return base;
    const pc = phaseColor(previewIx.phase);
    // Label each claimed peak with its ratio term (√2, √3, 2 …) so the hover
    // draws the eye to WHICH peaks the candidate is built on — the recolour
    // alone is too subtle on small markers.
    const labelByPeak = new Map<number, string>();
    for (const pr of previewIx.peaks) {
      const t = ratioTerm(previewIx.phase, pr.ratio_position);
      if (t) labelByPeak.set(pr.peak_id, t);
    }
    return {
      ...base,
      peaks: base.peaks.map((p) => {
        if (!previewClaim.has(p.id)) return p;
        const label = labelByPeak.get(p.id);
        return { ...p, color: pc, ...(label ? { label } : {}) };
      }),
    };
  }, [traceQ.data, peaks, phase, previewIx, previewClaim]);

  // Rings + their caption phases come from ONE adapter walk: `phases` names
  // only the phases that actually emitted a ring on this frame, so the caption
  // can never chip a hue that is not there (FO-RING single source; the
  // ring-less fully-landed custom index is the mainline case it guards).
  const { rings, phases: ringPhases } = useMemo(
    () => toDetectorRings(activeIndices, peaks),
    [activeIndices, peaks],
  );
  const { assigned: combAssigned, leftover: combLeftover } = useMemo(
    () => toCombSeries(activeIndices, peaks),
    [activeIndices, peaks],
  );

  // Detector image url for the active exposure (loupe's builder → cache-coherent).
  const detectorSrc = activeExposure
    ? buildExposureImageUrl(activeExposure)
    : null;

  const exposureLabel = activeExposure?.filename
    ? activeExposure.filename.replace(/\.[^.]+$/, "")
    : null;
  const sampleName =
    corpusSample?.name ?? "—";
  // FO-RESCORE2 F14: name the browser tab after the sample (was static
  // "Himalaya"). Raw name (null while loading) so the placeholder "—" never
  // leaks into the tab. Hooks-safe: above the not-found early return.
  useDocumentTitle(corpusSample?.name ?? null);
  // sanitizeDashes: upstream experiment/sample names may carry em dashes that
  // no source-level guard can catch (FO-SUBTITLE-EMDASH). Fold them out.
  const subtitle = sanitizeDashes(
    [corpusSample?.name, experimentQ.data?.name, exposureLabel]
      .filter(Boolean)
      .join(" · "),
  );

  // custom-index live preview
  const customMeta = CUSTOM_SYMS.find((s) => s.name === customSym) ?? CUSTOM_SYMS[0]!;
  const observedQs = peaks.filter((p) => !p.excluded).map((p) => p.q);
  const { previewSeries, fit } = customIndexPreview(
    customSym,
    Number(customParam),
    observedQs,
  );

  // Client-side validation of the lattice parameter, mirroring the trace q-add
  // field: an empty / non-finite / out-of-range value disables "Add to
  // assignment" so it never round-trips to a backend 400. (FO-PARAM-VALIDATION)
  const customParamNum = Number(customParam);
  const customParamValid =
    customParam.trim() !== "" &&
    Number.isFinite(customParamNum) &&
    customParamNum >= customMeta.min &&
    customParamNum <= customMeta.max;

  // FO-NAV-SKELETON: a sample switch ([ ]) clears activeExposureId (the
  // setActiveSample cascade) BEFORE the new sample's exposures round-trip and
  // useAutoPickExposure re-seeds it. In that window the trace/peaks queries are
  // disabled, so gating purely on `traceQ.isLoading` behind the
  // `activeExposureId !== undefined` guard read false and flashed empty panels
  // before the boneyard appeared. Treat the whole exposure-resolution window as
  // loading: we're still resolving WHICH exposure to show whenever there's no
  // active exposure AND the sample isn't confirmed to have zero usable ones
  // (that case is the honest "no exposures" empty state, not a skeleton). The
  // window ends when an exposure resolves (→ the trace/peaks check takes over).
  const { noUsable: noUsableExposure } = noUsableExposureState(exposuresQ.data);
  const resolvingExposure = activeExposureId === undefined && !noUsableExposure;
  const isLoading =
    corpusQ.isLoading ||
    resolvingExposure ||
    (activeExposureId !== undefined && (traceQ.isLoading || peaksQ.isLoading));

  // ── focus the scope when the skeleton reveals content (cold-load arrow nav) ──
  // The scope div attaches while the boneyard overlay is still visibility:hidden,
  // so el.focus() in the callback ref is silently dropped by Chromium. This effect
  // fires once per sample AFTER isLoading goes false (skeleton reveals the real
  // content), guaranteeing the scope is focusable before the first arrow keydown.
  useEffect(() => {
    if (isLoading) return;
    if (focusedSampleRef.current === activeSampleId) return;
    const el = scopeEl.current;
    if (el) { focusedSampleRef.current = activeSampleId ?? null; el.focus({ preventScroll: true }); }
  }, [isLoading, activeSampleId]);

  // FO-COMB-AXIS: the q-domain the comb shares with the trace. A manual zoom
  // (`xDomain`) wins; otherwise the trace auto-fits to its data extent, so the
  // comb mirrors that positive-q extent. Keeps the comb's axis ticks + reflection
  // positions in the same q-space the trace shows above it.
  const combSharedDomain = useMemo<[number, number] | null>(() => {
    if (xDomain) return xDomain;
    const qs = traceQ.data?.q;
    if (!qs || qs.length === 0) return null;
    let lo = Infinity;
    let hi = -Infinity;
    for (const q of qs) {
      if (q > 0) {
        if (q < lo) lo = q;
        if (q > hi) hi = q;
      }
    }
    return Number.isFinite(lo) && hi > lo ? [lo, hi] : null;
  }, [xDomain, traceQ.data]);

  // ── figure export ─────────────────────────────────────────────────────────────
  // The single-trace Focus figure renders through the SAME greenfield builder as
  // the series waterfall (buildCleanFigureSvg): build the one assignment row, then
  // render it in the clean export skin. WYSIWYG holds — it reads the live trace,
  // the active assignment's claimed peaks, the assigned phase, and the log/linear
  // toggle. (Was the legacy Observable Plot renderer, now fully retired.)
  const figureRow = useMemo(
    () =>
      buildFocusFigureRow({
        trace: traceQ.data ?? EMPTY_TRACE,
        peaks,
        activeIndices,
        phase,
        // No right-gutter row label on a single-trace figure: the title already
        // carries sample · exposure, so a gutter dup just clips against the
        // margin. (The gutter earns its keep only in the multi-row waterfall.)
        label: "",
      }),
    [traceQ.data, peaks, activeIndices, phase],
  );
  const figureTitle = [sampleName, exposureLabel].filter(Boolean).join(" · ");

  // The single-trace figure's KEY (legend): the trace's phase colour + every
  // active-assignment phase with its lattice a/d and, for cubics, κ. Coexisting
  // phases each get a segment (equal billing). The label is empty — the title
  // already carries sample · exposure, so the key reads as a phase/a/κ legend.
  const figureKey = useMemo<FigureTraceKey>(() => {
    const qUnits = experimentQ.data?.q_units;
    const seen = new Set<string>();
    const segments = [];
    for (const ix of activeIndices) {
      if (seen.has(ix.phase)) continue;
      seen.add(ix.phase);
      segments.push(phaseSegment(ix.phase, ix.lattice_d, qUnits));
    }
    const key: FigureTraceKey = { color: phaseHex(phase), label: "", segments };
    if (segments.length === 0) {
      key.note = assignment?.state === "form_factor" ? "form factor (no Bragg)" : "unindexed";
    }
    return key;
  }, [activeIndices, phase, experimentQ.data, assignment?.state]);

  // Descriptive, product-tagged stem (buildFilename slugifies it): e.g.
  // "himalaya-trace-jc042-frame-1-2026-06-13.svg".
  const filenameStem = `himalaya-trace-${sampleName} ${exposureLabel ?? ""}`.trim();
  const renderSvg = useCallback(
    () =>
      buildCleanFigureSvg({
        rows: [figureRow],
        traceKeys: [figureKey],
        title: figureTitle || "Trace",
        footer: experimentQ.data?.name ?? "",
        xType: scale === "log" ? "log" : "linear",
        qDomain: xDomain ?? waterfallQDomain([figureRow]),
        showPeakLabels: true,
      }),
    [figureRow, figureKey, figureTitle, experimentQ.data, scale, xDomain],
  );
  const fx = useFigureExport(renderSvg, filenameStem, "trace plot");

  // ── early states ─────────────────────────────────────────────────────────────
  if (routeStatus === "unknown" || (!corpusQ.isLoading && !corpusSample)) {
    return (
      <div className="px-8 py-7">
        <div data-testid="focus-not-found">
          <EmptyState
            as="h1"
            title="Sample not found"
            body="Nothing in the corpus matches this address."
            action={
              <Button variant="outline" onClick={() => navigate("/experiments")}>
                Back to the experiments
              </Button>
            }
          />
        </div>
      </div>
    );
  }

  // The honest empty state: the sample's exposures have RESOLVED to zero usable
  // ones (loaded, all rejected or none) — not the transient resolving window
  // where activeExposureId is momentarily undefined but acceptable exposures
  // exist (that window is the skeleton's, gated by `resolvingExposure` above).
  const noExposure = noUsableExposure;

  // ── custom-index helper ──────────────────────────────────────────────────────
  function commitCustom(): void {
    if (!customParamValid) return; // defense in depth; the Add button is disabled too
    commitCustomIndex.mutate(customSym, basisFor(customSym, Number(customParam)));
    setCustomOpen(false);
    // Consequential: a custom hypothesis was added to the call → visible toast.
    showToast(`${customSym} index added`, "success");
  }

  // ── rail building ────────────────────────────────────────────────────────────
  // Speculative-kind indices are the output of the custom-index builder (it
  // persists a speculative index + adds it to the assignment); they render as
  // ordinary toggle candidate rows alongside the auto candidates — the mockup's
  // candidate list is uniform/toggle-only (no per-row delete). The custom-index
  // modal is the hypothesis tool; the cart's remove (PhaseBlock onRemove) takes
  // a phase back out of the call.
  const speculatives = indices.filter((i) => i.kind === "speculative");

  // Single shared predicate: the cart's empty copy AND the candidate-list line
  // both branch on peaksEmpty, so the two texts can never contradict each
  // other (with zero peaks, "Every peak is unindexed" / "Check a candidate
  // below" would both be lies — the real next action is marking peaks).
  const peaksEmpty = peaks.length === 0;
  // FO-ALLEXCLUDED-CAPTION: peaks exist but every one is excluded → there is
  // nothing indexable, which is NOT the same as "candidates were tried and none
  // fit". Give it its own honest copy (and keep the candidate-list line below in
  // lockstep) so the cart never implies a failed search over an empty input.
  const allExcluded = !peaksEmpty && peaks.every((p) => p.excluded);
  const cartEmpty = peaksEmpty
    ? "No peaks marked. Find peaks on the trace to start indexing."
    : allExcluded
      ? "All peaks are excluded. Restore a peak, or add one, to index."
      : candidatePool.length === 0 && speculatives.length === 0
        ? "No phase assigned. No candidate fits these peaks. Try a custom index."
        : undefined; // candidates exist → AssignmentCart's default copy is correct

  function phaseMeta(ix: IndexEntry): string {
    const lattice =
      ix.lattice_d != null ? `a = ${ix.lattice_d.toFixed(0)} Å · ` : "";
    return `${lattice}${ix.peaks.length} reflections`;
  }

  const cart = (
    <AssignmentCart
      onCustomIndex={() => setCustomOpen(true)}
      {...(cartEmpty !== undefined ? { empty: cartEmpty } : {})}
    >
      {activeIndices.map((ix) => (
        <PhaseBlock
          key={ix.id}
          phase={ix.phase}
          score={ix.score ?? 0}
          meta={phaseMeta(ix)}
          series={seriesRatio(ix.phase, ix.peaks.map((p) => p.ratio_position)) || undefined}
          onRemove={() => {
            removeAssignmentPhase.mutate(ix.id);
            announce(`${ix.phase} removed from the call`);
          }}
        />
      ))}
    </AssignmentCart>
  );

  function candidateRow(ix: IndexEntry): JSX.Element {
    const selected = memberIds.has(ix.id);
    return (
      // Placement-only wrapper: hover → transient hoverPreviewId; click → set
      // cursor + sticky preview (keyboard position follows pointer).
      // Intentional dual-fire: the wrapper onClick sets the cursor position while
      // CandidateRow's onClick (onToggle) toggles assignment. CandidateRow must
      // NOT call e.stopPropagation() — doing so would silently break cursor-setting.
      <div
        key={ix.id}
        onMouseEnter={() => setHoverPreviewId(ix.id)}
        onMouseLeave={() => setHoverPreviewId(undefined)}
        onClick={() => { candidateCursor.setCursor(ix.id); setPreviewWasExplicit(true); }}
      >
        <CandidateRow
          phase={ix.phase}
          score={ix.score}
          why={`explains ${ix.peaks.length} peaks${selected ? " · in the call" : ""}`}
          selected={selected}
          previewed={ix.id === previewIndexId}
          {...(ix.bonnet?.consistent ? { bonnet: true } : {})}
          onToggle={() => {
            if (selected) {
              removeAssignmentPhase.mutate(ix.id);
              announce(`${ix.phase} removed from the call`);
            } else {
              addAssignmentPhase.mutate(ix.id);
              announce(`${ix.phase} added to the call`);
            }
          }}
        />
      </div>
    );
  }

  // The form-factor declaration ("no Bragg peaks to index") is the honest
  // alternative to indexing — it only makes sense while nothing is in the call,
  // so it is hidden the moment a phase is confirmed (it would otherwise clear
  // that phase on the next click). It rides the foot of the candidate list so
  // declaring it is the same one-check gesture as picking a phase, but a
  // hairline sets it apart because it is not a phase. Checking a phase candidate
  // flips the assignment back to `indexed` (the add mutator), so the two are
  // mutually exclusive without any extra wiring here.
  const isFormFactor = assignment?.state === "form_factor";
  const candidates = (
    <>
      <CandidateList>
        {candidatePool.length === 0 ? (
          <HintText>
            {peaksEmpty
              ? "Candidates appear once peaks are marked."
              : allExcluded
                ? "Candidates appear once a peak is restored."
                : "No candidate indexings."}
          </HintText>
        ) : (
          candidatePool.map((ix) => candidateRow(ix))
        )}
        {speculatives.map((ix) => candidateRow(ix))}
      </CandidateList>
      {activeIndices.length === 0 && (
        <div className="mt-3 border-t border-hair pt-3">
          <FormFactorRow
            selected={isFormFactor}
            onToggle={() => {
              setAssignmentState.mutate(isFormFactor ? "null" : "form_factor");
              announce(isFormFactor ? "Form factor cleared" : "Marked as form factor");
            }}
          />
        </div>
      )}
    </>
  );

  // ── modals (mounted regardless of layout branch) ─────────────────────────────
  const modals = (
    <>
      <CustomIndexModal
        open={customOpen}
        symmetries={CUSTOM_SYMS.map((s) => s.name)}
        symmetry={customSym}
        onSymmetryChange={setCustomSym}
        paramName={customMeta.paramName}
        paramMin={customMeta.min}
        paramMax={customMeta.max}
        {...(customMeta.step !== undefined ? { paramStep: customMeta.step } : {})}
        unit={customMeta.unit}
        paramValue={customParam}
        onParamChange={setCustomParam}
        previewSeries={previewSeries}
        observed={observedQs}
        fit={fit}
        addDisabled={!customParamValid}
        paramInvalid={customParam.trim() !== "" && !customParamValid}
        onAdd={commitCustom}
        onCancel={() => setCustomOpen(false)}
        onClose={() => setCustomOpen(false)}
      />
    </>
  );

  // ── layout ───────────────────────────────────────────────────────────────────
  // Full-bleed per the focus-plot mockup: a [work 1fr · rail 350px] grid with
  // NO outer max-width — the work column centres its own content at 1180px and
  // the rail (bg-paper-sunk + left hairline) pins to the right edge.
  return (
    <>
      <Skeleton
        name="focus"
        className="block"
        loading={isLoading}
        stagger={50}
        transition={200}
        fixture={FOCUS_FIXTURE}
        fallback={
          <div
            data-testid="focus-skeleton"
            className="p-8 text-sm italic text-ink-soft"
          >
            Loading workspace…
          </div>
        }
      >
        <div
          ref={scopeRef}
          tabIndex={-1}
          data-interaction-scope
          data-testid="focus-workspace"
          className={FOCUS_PAGE_GRID}
        >
          {/* work column — full-bleed; inner content capped at 1180px (mockup .work / .work-inner) */}
          <div className="min-w-0 px-8 pt-7 pb-13">
           <div className="mx-auto flex min-w-0 max-w-[1180px] flex-col gap-[18px]">
            <TracePlate
              kicker="Integration"
              title={sampleName}
              subtitle={subtitle}
              trace={traceModel}
              scale={scale}
              onScaleChange={setScale}
              // The Focus curve stays neutral gray (an assignment hue on the
              // trace is arbitrary under coexistence); only the peak markers and
              // their hover labels carry phase colour. Headroom lifts the
              // ceiling so the tallest peak's marker clears the curve. Ratio
              // labels show only while previewing a candidate.
              neutralLine
              yHeadroom={0.4}
              showPeakLabels={previewClaim !== undefined}
              addPeakArmed={addArmed}
              onToggleAddPeak={() => setAddArmed((v) => !v)}
              onAutoFit={() => setXDomain(null)}
              xDomain={xDomain}
              {...(hoveredQ !== undefined ? { hoveredQ } : {})}
              onHoverQ={setHoveredQ}
              {...(previewClaim !== undefined ? { highlightPeakIds: previewClaim } : {})}
              {...(focusRequest !== undefined ? { focusRequest } : {})}
              interaction={{
                onXDomain: setXDomain,
                onAddPeak: (q) => {
                  addPeak.mutate(q);
                  announce("Peak added");
                },
                onClickPeak: (id) => {
                  // Provenance decides the verb (peakClickAction): an auto peak
                  // is the indexer's — a click toggles its excluded flag (disable
                  // / restore), which restyles the mark in place WITHOUT
                  // unmounting it, so focus stays put. A manual peak is the
                  // user's — a click removes it, unmounting the activated mark, so
                  // re-anchor focus to the surviving q-neighbour (computed from
                  // the pre-removal list; "+ Peak" button via fallback if none).
                  const action = peakClickAction(peaks, id);
                  if (!action) return;
                  if (action.kind === "toggle-exclude") {
                    setPeakExcluded.mutate({ peakId: id, excluded: action.excluded });
                    announce(action.excluded ? "Auto peak disabled" : "Auto peak restored");
                  } else {
                    requestPeakFocus(nextFocusPeakId(peaks, id));
                    removePeak.mutate(id);
                    announce("Peak removed");
                  }
                },
                onReset: () => setXDomain(null),
              }}
              actions={
                <ExportButton
                  {...fx}
                  ariaContext="trace plot"
                  // Export gates on trace emptiness only: a peakless trace is
                  // still an exportable WYSIWYG figure (BU-EXPORT precedent).
                  disabled={!traceQ.data || traceQ.data.q.length === 0}
                />
              }
            />
            <div className={FOCUS_SPLIT_GRID}>
              <DetectorPanel
                src={detectorSrc}
                rings={rings}
                ringPhases={ringPhases}
                orient={detOrient}
                onRawSize={handleRawSize}
                onOrient={setDetOrient}
                {...(calibration ? { calibration } : {})}
                {...(imageAspect !== undefined ? { imageAspect } : {})}
                {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                onHoverQ={setHoveredQ}
                tools={
                  exposures.length > 0 ? (
                    // FO-DETSTRIP-WALL: render at "sm" (62px), not "xs" (30px) —
                    // at 30px the detector content was an indistinguishable dark
                    // square, the frame-number label ate a third of the tile, and
                    // the selected ring read as a hairline. 62px shows real
                    // per-frame content + a legible frame label + a clear ring.
                    <ThumbnailGallery
                      exposures={toGalleryExposures(exposures)}
                      {...(activeExposureId !== undefined ? { selectedId: activeExposureId } : {})}
                      onSelect={setActiveExposure}
                      size="sm"
                    />
                  ) : undefined
                }
              />
              <div className="hidden lg:flex min-h-0 flex-col">
                <CombsPanel
                  assigned={combAssigned}
                  leftover={combLeftover}
                  view={combView}
                  onViewChange={setCombView}
                  {...(combSharedDomain ? { qDomain: combSharedDomain } : {})}
                  {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                  onHoverQ={setHoveredQ}
                />
              </div>
            </div>
           </div>
          </div>

          {/* rail column = AssignmentRail itself (it owns the rail shell:
              bg-paper-sunk + left hairline + padding + scroll). The empty state
              mirrors that shell so the column looks identical with no exposures. */}
          {noExposure ? (
            <div className="flex min-h-0 flex-col border-l border-hair bg-paper-sunk p-5 text-sm text-ink-soft">
              This sample has no exposures.
            </div>
          ) : (
            <AssignmentRail
              assignmentCount={activeIndices.length || undefined}
              assignment={cart}
              candidates={candidates}
              candidatesNote="A sample can be multiphasic, so check every phase that fits."
            />
          )}
        </div>
      </Skeleton>

      {modals}
    </>
  );
}
