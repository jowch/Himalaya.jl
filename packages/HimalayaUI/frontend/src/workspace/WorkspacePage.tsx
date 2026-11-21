import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import type { Exposure, IndexEntry } from "../api";
import { ApiError } from "../api";
import {
  useAddAssignmentPhase,
  useAddPeak,
  useCommitCustomIndex,
  useAssignment,
  useBridgeSamples,
  useBridges,
  useExperimentCorpusSamples,
  useExperiment,
  useExposures,
  useIndices,
  useLabels,
  useAddLabel,
  useRemoveLabel,
  useProblems,
  useDeleteExperiment,
  useLabelColors,
  useSetLabelColor,
  usePeaks,
  useFrameAssignments,
  useRemoveAssignmentPhase,
  useRemoveFramePhase,
  useSampleCall,
  useSetSampleCallState,
  useRemovePeak,
  useSetExposureStatus,
  useSetPeakExcluded,
  useTrace,
  useUnbridgeSample,
  useTriggerScan,
  useCheckForNewFilesOnOpen,
} from "../queries";
import { useScanState } from "../hooks/useScanState";
import { useMediaQuery } from "./useMediaQuery";
import { isIngestStage, stageLabel } from "../lib/ingestStages";
import { bridgeCounts, bridgedWith } from "../lib/bridges";
import { bufferRole, buffersFirst, displayLabel, labelsFor, sampleMatchesFilter, labelVocabulary } from "../lib/labels";
import { representativeExposure } from "../lib/exposure";
import { phaseColor } from "../phases";
import { useAppState } from "../state";
import { MAX_CANDIDATES, firstUnreviewed, i0BelowZero, traceIsDark, imageUrl, nextUnreviewed, peakClick, rankCandidates, sampleStatus, stepExposure } from "./model";
import { samePeak, type PeakRef } from "../lib/peakRef";
import { addCompare, clampDomain, compareScale, scaledTrace, subtractTrace } from "./compare";
import type { CompareScale } from "./compare";
import { useCompareTraces } from "./useCompareTraces";
import { usePrefetchNeighbours } from "./prefetch";
import { CompareBar } from "./ui/CompareBar";
import { ComparePicker } from "./ui/ComparePicker";
import { SubtractionStrip } from "./ui/SubtractionStrip";
import { PlotCard } from "./ui/PlotCard";
import { HelpOverlay, type HelpSection } from "./ui/HelpOverlay";
import { DetectorOverlay } from "./ui/DetectorOverlay";
import { Palette } from "./ui/palette/Palette";
import { workspacePalette } from "./palette";
import { useWorkspaceTheme } from "./theme";
import { PLOT_VIEWS, usePlotView } from "./plotView";
import { positiveExtent } from "../print/plot/projection";
import { CustomIndexPanel, NUDGE_STEP, NUDGE_STEP_LARGE } from "./ui/CustomIndexPanel";
import { IndexingPlot, type IndexSeries } from "./ui/IndexingPlot";
import { clampLattice, draftCommit, draftPrediction, indexFit, initialDraft, storedFit, type CustomDraft } from "./indexing";
import { latticeBounds, latticeForFirstOrderOnPeak } from "../lib/customIndex";

import { WorkspaceHeader } from "./ui/WorkspaceHeader";
import { Note } from "./ui/Note";
import { WorkspaceLayout } from "./ui/WorkspaceLayout";
import { SampleList } from "./ui/SampleList";
import { ExperimentNotice } from "./ui/ExperimentNotice";
import { ProblemList } from "./ui/ProblemList";
import { RemoveExperimentDialog } from "./ui/RemoveExperimentDialog";
import { COMPARE_COLORS, WorkspacePlot, type CombLine, type CompareLine } from "./ui/WorkspacePlot";
import { ExposureStrip } from "./ui/ExposureStrip";
import { CandidateList } from "./ui/CandidateList";
import { MetaList } from "./ui/MetaList";
import { MainSection, SampleHeading, Section } from "./ui/Section";
import { CustomIndexButton } from "./ui/CustomIndexButton";
import { LabelRow } from "./ui/LabelRow";
import { SYMS } from "../lib/customIndex";

function helpSections(mod: string): HelpSection[] {
  return [
    {
      title: "Move",
      rows: [
        { keys: ["↑", "↓"], what: "Previous or next sample in the list" },
        { keys: ["Z"], what: "Next unreviewed sample" },
        { keys: ["/"], what: "Filter samples by name or label" },
      ],
    },
    {
      title: "Exposures",
      rows: [
        { keys: ["←", "→"], what: "Previous or next exposure" },
        { keys: ["X"], what: "Drop or restore the exposure" },
        { keys: ["Space"], what: "Enlarge its detector image; ← → and X work inside it" },
      ],
    },
    {
      title: "Call",
      rows: [
        { keys: ["1–9"], what: "Tick or untick a phase candidate" },
        { keys: ["0"], what: "Mark no phase (form factor), or clear it" },
        { keys: ["I"], what: "Custom index: pick a symmetry and drag the lattice; Enter adds, Esc cancels" },
        { keys: ["L"], what: "Label the sample, e.g. Buffer; Enter adds, Esc closes" },
      ],
    },
    {
      title: "Plot",
      rows: [
        { gesture: "Click", what: "Turn a peak off or on, remove one you added, or add one" },
        { gesture: "Scroll", what: "Zoom q" },
        { keys: [`${mod} 0`], what: "Reset the zoom" },
        { keys: ["C"], what: "Compare with another sample" },
      ],
    },
    {
      title: "Anywhere",
      rows: [
        { keys: [mod, "K"], what: "Commands and go to a sample" },
        { keys: ["?"], what: "This list" },
      ],
    },
  ];
}

const MOD_KEY = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.userAgent) ? "⌘" : "Ctrl";
const HELP = helpSections(MOD_KEY);

function latticeParam(phase: string): string {
  return SYMS[phase]?.param ?? "a";
}

function exposureRows(e: Exposure): [string, string][] {
  const rows: [string, string][] = [["File", e.filename ?? "—"]];
  rows.push(["Kind", e.kind === "averaged" ? "Averaged" : e.kind === "background_subtracted" ? "Subtracted" : "Single"]);
  const i0 = e.incident_monitor;
  const i2 = e.transmitted_counts;
  rows.push(["Acquired", e.timestamp ? e.timestamp.replace("T", " ").slice(0, 19) : "—"]);
  rows.push(["Exposure", e.exposure_time != null ? `${e.exposure_time.toFixed(1)} s` : "—"]);
  rows.push(["Flight tube", e.pipe_length_m != null ? `${Math.round(e.pipe_length_m * 1000)} mm` : "—"]);
  rows.push(["Energy", e.beam_energy_ev != null ? `${(e.beam_energy_ev / 1000).toFixed(3)} keV` : "—"]);
  rows.push(["Detector", e.detector || "—"]);
  rows.push(["I₀", i0 != null ? i0.toLocaleString("en-US", { maximumSignificantDigits: 4 }) : "—"]);
  rows.push(["I signal", i2 != null ? Math.round(i2).toLocaleString("en-US") : "—"]);
  return rows;
}

const SHORT_WINDOW = "(max-height: 860px)";

function isTyping(target: EventTarget | null): boolean {
  return target instanceof HTMLElement && (target.isContentEditable || ["INPUT", "TEXTAREA", "SELECT"].includes(target.tagName));
}

function inOverlay(target: EventTarget | null): boolean {
  return target instanceof HTMLElement && target.closest('[role="dialog"],[data-keys-trap]') !== null;
}

export function WorkspacePage(): JSX.Element {
  const { id, sampleId: sampleParam } = useParams<{ id: string; sampleId?: string }>();
  const experimentId = Number(id);
  const navigate = useNavigate();
  const [search] = useSearchParams();
  const linkedParam = Number(search.get("exposure"));
  const linkedExposure = Number.isInteger(linkedParam) && linkedParam > 0 ? linkedParam : undefined;
  const username = useAppState((s) => s.username);
  const clearUsername = useAppState((s) => s.clearUsername);
  const experiment = useExperiment(experimentId);
  const scan = useScanState(experimentId);
  const inFlight = useAppState((s) => s.ingestInFlight?.[experimentId]);
  const triggerScan = useTriggerScan(experimentId);
  useCheckForNewFilesOnOpen(experimentId);
  const problems = useProblems(experimentId);
  const [problemsOpen, setProblemsOpen] = useState(false);
  const [removeOpen, setRemoveOpen] = useState(false);
  const removeExperiment = useDeleteExperiment(experimentId);
  const problemRows = problems.data?.exposures ?? [];
  const waitingForTrace = problemRows.filter((p) => p.problem === "no_trace").length;
  const failedAnalyses = problemRows.length - waitingForTrace;
  const ingesting = scan.scanning || scan.rescanning;
  const count = (n: number): string => n.toLocaleString("en-US");
  const ingestProgress = !ingesting
    ? undefined
    : inFlight && inFlight.total > 0
      ? `${isIngestStage(inFlight.stage) ? stageLabel(inFlight.stage) : "Ingesting"} ${count(inFlight.processed)} of ${count(inFlight.total)}`
      : scan.scanning ? "Reading files…" : "Checking for new files…";
  const ingestFailed = scan.failed || triggerScan.isError;
  const corpus = useExperimentCorpusSamples(experimentId);
  const [filter, setFilter] = useState("");

  const samples = useMemo(
    () => (corpus.data ?? []).filter((s) => s.experiment_id === experimentId),
    [corpus.data, experimentId],
  );
  const labels = useLabels(experimentId);
  const labelList = useMemo(() => labels.data ?? [], [labels.data]);
  const vocabulary = useMemo(() => labelVocabulary(labelList), [labelList]);
  const addLabel = useAddLabel(experimentId);
  const removeLabel = useRemoveLabel(experimentId);
  const [labeling, setLabeling] = useState(false);
  const labelColors = useLabelColors();
  const setLabelColor = useSetLabelColor();
  const [coloringLabel, setColoringLabel] = useState<string | undefined>(undefined);
  const visible = useMemo(
    () => samples.filter((s) => sampleMatchesFilter(s.name, labelsFor(labelList, s.id), filter)),
    [samples, labelList, filter],
  );

  const requested = sampleParam === undefined ? undefined : Number(sampleParam);
  const sample = requested === undefined ? firstUnreviewed(samples) : samples.find((s) => s.id === requested);
  const sampleId = sample?.id;
  const goSample = useCallback(
    (sid: number) => navigate(`/experiments/${experimentId}/workspace/${sid}`, { replace: true }),
    [navigate, experimentId],
  );
  useEffect(() => {
    if (requested !== undefined || sampleId === undefined) return;
    const qs = search.toString();
    navigate(`/experiments/${experimentId}/workspace/${sampleId}${qs ? `?${qs}` : ""}`, { replace: true });
  }, [requested, sampleId, search, navigate, experimentId]);

  const exposures = useExposures(sampleId);
  const [exposureId, setExposureId] = useState<number | undefined>(linkedExposure);
  useEffect(() => setExposureId(linkedExposure), [sampleId, linkedExposure]);
  const exposureList: Exposure[] = useMemo(() => exposures.data ?? [], [exposures.data]);
  const exposure = exposureList.find((e) => e.id === exposureId) ?? representativeExposure(exposureList);

  const trace = useTrace(exposure?.id);
  const noTrace = trace.error instanceof ApiError && trace.error.status === 404;
  const i0Hint = i0BelowZero(exposure);
  const darkNotice = traceIsDark(trace.data)
    ? `No counts on the detector: the beam was likely off or the shutter closed.${i0Hint ? " I₀ also reads below zero." : ""}`
    : undefined;
  const noTraceMessage = noTrace
    ? `No integrated trace for this frame.${i0Hint ? " I₀ reads below zero, so the beam was likely off." : ""}`
    : undefined;
  const peaks = usePeaks(exposure?.id);
  const indices = useIndices(exposure?.id);
  const assignment = useAssignment(exposure?.id);
  usePrefetchNeighbours({
    samples: visible,
    currentId: sample?.id,
    exposures: exposureList,
    shownId: exposure?.id,
    settled: [trace, peaks, indices, assignment].every((q) => q.fetchStatus === "idle"),
  });
  const allCandidates = useMemo(() => rankCandidates(indices.data, Infinity), [indices.data]);
  const candidates = useMemo(() => allCandidates.slice(0, MAX_CANDIDATES), [allCandidates]);
  const sampleCall = useSampleCall(sampleId);
  const otherFrames = useMemo(
    () => [...new Set((sampleCall.data?.members ?? []).map((m) => m.exposure_id))].filter((id) => id !== exposure?.id),
    [sampleCall.data, exposure?.id],
  );
  const otherAssignments = useFrameAssignments(otherFrames);
  const otherMembers = useMemo(() => {
    const live = new Map(otherFrames.map((id, i) => [id, otherAssignments[i]?.data?.members]));
    return (sampleCall.data?.members ?? []).filter(
      (m) => m.exposure_id !== exposure?.id && (live.get(m.exposure_id) ?? [m.index_id]).includes(m.index_id),
    );
  }, [sampleCall.data, otherFrames, otherAssignments, exposure?.id]);
  const assigned = useMemo(() => {
    const here = new Set(assignment.data?.members ?? []);
    const phases = new Set(otherMembers.map((m) => m.phase));
    for (const c of allCandidates) if (phases.has(c.phase)) here.add(c.id);
    return here;
  }, [assignment.data, otherMembers, allCandidates]);
  const [hoverId, setHoverId] = useState<number | undefined>(undefined);

  const addPhase = useAddAssignmentPhase(exposure?.id ?? 0);
  const removePhase = useRemoveAssignmentPhase(exposure?.id ?? 0);
  const removeFramePhase = useRemoveFramePhase();
  const setCallState = useSetSampleCallState(sampleId ?? 0);
  const isFormFactor = assignment.data?.state === "form_factor";
  const setStatus = useSetExposureStatus(sampleId ?? 0);

  const exposureDropped = exposure?.status === "rejected";
  const canEditPeaks = Boolean(username && exposure && !exposureDropped);
  const addPeak = useAddPeak(exposure?.id ?? 0);
  const removePeak = useRemovePeak(exposure?.id ?? 0);
  const setExcluded = useSetPeakExcluded(exposure?.id ?? 0);

  const [xDomain, setXDomain] = useState<[number, number] | null>(null);
  const [view, setView] = usePlotView();
  const [refs, setRefs] = useState<number[]>([]);
  const [subChoice, setSubRef] = useState<number | undefined>(undefined);
  const [refsOn, setRefsOn] = useState(true);
  const [picking, setPicking] = useState(false);
  const [helpOpen, setHelpOpen] = useState(false);
  const [bigOpen, setBigOpen] = useState(false);
  const [paletteOpen, setPaletteOpen] = useState(false);
  const theme = useWorkspaceTheme();
  const [draft, setDraft] = useState<CustomDraft | null>(null);
  useEffect(() => setDraft(null), [exposure?.id]);
  const commitCustom = useCommitCustomIndex(exposure?.id ?? 0);
  const bridges = useBridges(experimentId);
  const bridgeList = useMemo(() => bridges.data ?? [], [bridges.data]);
  const bridged = useMemo(() => (sampleId === undefined ? [] : bridgedWith(bridgeList, sampleId)), [bridgeList, sampleId]);
  const counts = useMemo(() => bridgeCounts(bridgeList), [bridgeList]);
  const bridgeSamples = useBridgeSamples(experimentId);
  const unbridgeSample = useUnbridgeSample(experimentId);
  const shownRefs = useMemo(() => [...bridged, ...refs.filter((r) => !bridged.includes(r))], [bridged, refs]);
  const subRef = subChoice !== undefined && shownRefs.includes(subChoice) ? subChoice : undefined;
  const compareTraces = useCompareTraces(shownRefs);

  const onClickPeak = useCallback(
    (ref: PeakRef) => {
      const action = peakClick(peaks.data ?? [], ref);
      if (!action) return;
      if (action.kind === "remove") removePeak.mutate(ref.id);
      else setExcluded.mutate({ peakId: ref.id, excluded: action.excluded });
    },
    [peaks.data, removePeak, setExcluded],
  );

  const removeCompare = useCallback((rid: number) => {
    setRefs((r) => r.filter((x) => x !== rid));
    setSubRef((cur) => (cur === rid ? undefined : cur));
  }, []);

  const toggleCandidate = useCallback(
    (ix: IndexEntry) => {
      if (!exposure) return;
      if (assignment.data?.members.includes(ix.id)) {
        removePhase.mutate(ix.id);
        return;
      }
      const elsewhere = otherMembers.filter((m) => m.phase === ix.phase);
      if (elsewhere.length > 0) {
        for (const m of elsewhere) removeFramePhase.mutate(m.exposure_id, m.index_id);
        return;
      }
      if (exposure.status === "rejected") return;
      addPhase.mutate(ix.id);
    },
    [exposure, assignment.data, otherMembers, addPhase, removePhase, removeFramePhase],
  );

  const toggleFormFactor = useCallback(() => {
    if (!exposure || exposure.status === "rejected") return;
    setCallState.mutate(isFormFactor ? "null" : "form_factor");
  }, [exposure, isFormFactor, setCallState]);

  const toggleDrop = useCallback(
    (e: Exposure) => {
      setExposureId(e.id);
      setStatus.mutate({ exposureId: e.id, status: e.status === "rejected" ? null : "rejected" });
    },
    [setStatus],
  );

  const draftQs = useMemo(() => (draft ? draftPrediction(draft) : []), [draft]);
  const combs: CombLine[] = useMemo(() => {
    if (draft) return [{ id: -1, qs: draftQs, color: phaseColor(draft.phase) }];
    return candidates.filter((c) => c.id === hoverId).map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase) }));
  }, [candidates, hoverId, draft, draftQs]);

  const sampleNames = useMemo(() => new Map(samples.map((s) => [s.id, s.name])), [samples]);
  const activeRefs = shownRefs.filter((r) => r !== sampleId);
  const mainTrace = trace.data;
  const peakQs = useMemo(() => (peaks.data ?? []).map((p) => p.q), [peaks.data]);
  const scales = useMemo(() => {
    const out = new Map<number, CompareScale>();
    if (!mainTrace) return out;
    compareTraces.forEach((t, rid) => out.set(rid, compareScale(mainTrace, t, peakQs)));
    return out;
  }, [mainTrace, compareTraces, peakQs]);
  const compares: CompareLine[] = useMemo(() => {
    if (!mainTrace || !refsOn) return [];
    return shownRefs.flatMap((rid, i) => {
      const t = compareTraces.get(rid);
      const scale = scales.get(rid)?.scale;
      if (rid === sampleId || !t || scale === undefined) return [];
      const sc = scaledTrace(t, scale);
      return [{ id: rid, name: sampleNames.get(rid) ?? `Sample ${rid}`, q: sc.q, I: sc.I, color: COMPARE_COLORS[i % COMPARE_COLORS.length]!, emphasis: rid === subRef }];
    });
  }, [mainTrace, refsOn, shownRefs, compareTraces, scales, sampleId, sampleNames, subRef]);
  const subTrace = subRef !== undefined && subRef !== sampleId && refsOn ? compareTraces.get(subRef) : undefined;
  const subScale = subRef !== undefined ? scales.get(subRef)?.scale : undefined;
  const difference = useMemo(
    () => (mainTrace && subTrace && subScale !== undefined ? subtractTrace(mainTrace, subTrace, subScale) : null),
    [mainTrace, subTrace, subScale],
  );
  const xExtent = useMemo<[number, number] | null>(
    () => (mainTrace && mainTrace.q.length > 1 ? positiveExtent(mainTrace.q) : null),
    [mainTrace],
  );
  const shownDomain = clampDomain(xDomain, xExtent);
  const pickable = useMemo(
    () => buffersFirst(samples.filter((s) => s.id !== sampleId && !shownRefs.includes(s.id)), labelList),
    [samples, sampleId, shownRefs, labelList],
  );
  const bufferNotes = useMemo(() => {
    const notes = new Map<number, string>();
    for (const s of samples) {
      const role = bufferRole(labelList, s.id);
      if (role) notes.set(s.id, role === "confirmed" ? "Buffer" : "Suggested buffer");
    }
    return notes;
  }, [samples, labelList]);
  const bridgePickable = useMemo(() => samples.filter((s) => s.id !== sampleId && !bridged.includes(s.id)), [samples, sampleId, bridged]);
  useEffect(() => {
    setLabeling(false);
    setColoringLabel(undefined);
  }, [sampleId]);
  const sampleLabels = useMemo(() => (sampleId === undefined ? [] : labelsFor(labelList, sampleId)), [labelList, sampleId]);
  const openLabeling = useCallback(() => {
    if (!username || sampleId === undefined) return;
    setColoringLabel(undefined);
    setLabeling(true);
  }, [username, sampleId]);
  const leaveBridges = useCallback(() => {
    if (username && sampleId !== undefined && bridged.length > 0) unbridgeSample.mutate({ sampleId });
  }, [username, sampleId, bridged, unbridgeSample]);
  const onRemoveChip = useCallback(
    (rid: number) => {
      if (!bridged.includes(rid)) {
        removeCompare(rid);
        return;
      }
      if (username) unbridgeSample.mutate({ sampleId: rid });
    },
    [bridged, removeCompare, username, unbridgeSample],
  );

  const activePeakQs = useMemo(() => (peaks.data ?? []).filter((p) => !p.excluded).map((p) => p.q), [peaks.data]);
  const canCall = Boolean(username && exposure && !exposureDropped);
  const draftBounds = draft ? latticeBounds(draft.phase, xExtent) : null;
  const shortWindow = useMediaQuery(SHORT_WINDOW);
  const [openPanel, setOpenPanel] = useState<"call" | "exposure">("call");
  const callOpen = !shortWindow || openPanel === "call" || draft !== null;
  const exposureOpen = !shortWindow || (openPanel === "exposure" && draft === null);
  const openCustom = useCallback(() => {
    if (!canCall) return;
    setDraft((d) => d ?? initialDraft(candidates[0]));
  }, [canCall, candidates]);
  const setDraftClamped = useCallback(
    (d: CustomDraft) => setDraft({ phase: d.phase, lattice: clampLattice(d.phase, d.lattice, xExtent) }),
    [xExtent],
  );
  const addCustom = useCallback(() => {
    if (!draft || !exposure || !canCall) return;
    const c = draftCommit(draft);
    if (c.ratios.length === 0) return;
    commitCustom.mutate(c.phase, c.basis, c.ratios);
    setDraft(null);
  }, [draft, exposure, canCall, commitCustom]);
  const snapDraft = useCallback(
    (q: number) => {
      if (!draft) return;
      let target = q;
      let bestD = Infinity;
      for (const pq of activePeakQs) {
        const d = Math.abs(pq - q);
        if (d < bestD) {
          bestD = d;
          target = pq;
        }
      }
      setDraftClamped({ ...draft, lattice: latticeForFirstOrderOnPeak(draft.phase, target) });
    },
    [draft, activePeakQs, setDraftClamped],
  );
  const snapDraftToPeak = useCallback(
    (ref: PeakRef) => {
      const p = (peaks.data ?? []).find((pk) => samePeak(pk, ref));
      if (p) snapDraft(p.q);
    },
    [peaks.data, snapDraft],
  );
  const draftFit = useMemo(() => (draft ? indexFit(draftQs, activePeakQs, xExtent) : null), [draft, draftQs, activePeakQs, xExtent]);

  const indexSeries: IndexSeries[] = useMemo(() => {
    if (!xExtent) return [];
    if (draft && draftFit) {
      return [{ key: "draft", phase: draft.phase, color: phaseColor(draft.phase), fit: draftFit, lattice: draft.lattice, param: latticeParam(draft.phase), r2: draftFit.r2 }];
    }
    const hovered = candidates.find((c) => c.id === hoverId);
    const shown = hovered ? [hovered] : candidates.filter((c) => assigned.has(c.id));
    return shown.map((c) => {
      const fit = storedFit(c, xExtent);
      return { key: String(c.id), phase: c.phase, color: phaseColor(c.phase), fit, lattice: c.lattice_d, param: latticeParam(c.phase), r2: c.r_squared ?? fit.r2 };
    });
  }, [draft, draftFit, xExtent, candidates, hoverId, assigned]);

  useEffect(() => {
    const onKey = (ev: KeyboardEvent): void => {
      if ((ev.metaKey || ev.ctrlKey) && !ev.altKey && !ev.shiftKey && ev.key === "0" && !isTyping(ev.target) && !inOverlay(ev.target)) {
        ev.preventDefault();
        setXDomain(null);
        return;
      }
      if ((ev.metaKey || ev.ctrlKey) && !ev.altKey && !ev.shiftKey && ev.key.toLowerCase() === "k") {
        ev.preventDefault();
        setHelpOpen(false);
        setPicking(false);
        setLabeling(false);
        setColoringLabel(undefined);
        setPaletteOpen((o) => !o);
        return;
      }
      if (ev.defaultPrevented || ev.metaKey || ev.ctrlKey || ev.altKey) return;
      if (inOverlay(ev.target)) return;
      if (ev.key === "/" && !isTyping(ev.target)) {
        ev.preventDefault();
        document.getElementById("workspace-filter")?.focus();
        return;
      }
      if (isTyping(ev.target)) {
        if (ev.key === "Escape") (ev.target as HTMLElement).blur();
        return;
      }
      if (draft) {
        if (ev.key === "Escape") {
          ev.preventDefault();
          setDraft(null);
        } else if (ev.key === "Enter") {
          if (ev.target instanceof HTMLElement && ev.target.closest("button,a")) return;
          ev.preventDefault();
          addCustom();
        } else if (ev.key === "ArrowLeft" || ev.key === "ArrowRight") {
          ev.preventDefault();
          const step = (ev.key === "ArrowRight" ? 1 : -1) * (ev.shiftKey ? NUDGE_STEP_LARGE : NUDGE_STEP);
          setDraftClamped({ ...draft, lattice: draft.lattice + step });
        } else if (ev.key === "?") {
          ev.preventDefault();
          setHelpOpen(true);
        }
        return;
      }
      const idx = visible.findIndex((s) => s.id === sampleId);
      const step = (d: 1 | -1): void => {
        const next = visible[Math.min(visible.length - 1, Math.max(0, (idx < 0 ? 0 : idx + d)))];
        if (next) goSample(next.id);
      };
      let handled = true;
      switch (ev.key) {
        case "ArrowDown":
          step(1);
          break;
        case "ArrowUp":
          step(-1);
          break;
        case "ArrowRight":
        case "ArrowLeft": {
          const eid = stepExposure(exposureList, exposure?.id, ev.key === "ArrowRight" ? 1 : -1);
          if (eid !== undefined) setExposureId(eid);
          break;
        }
        case "c":
        case "C":
          setPicking(true);
          break;
        case "i":
        case "I":
          openCustom();
          break;
        case "l":
        case "L":
          if (username && sample) openLabeling();
          else handled = false;
          break;
        case "?":
          setHelpOpen(true);
          break;
        case "z":
        case "Z": {
          const nid = nextUnreviewed(visible, sampleId, 1);
          if (nid !== undefined) goSample(nid);
          break;
        }
        case " ":
          if (ev.repeat) break;
          if (exposure && imageUrl(exposure, false)) setBigOpen(true);
          else handled = false;
          break;
        case "x":
        case "X":
          if (exposure && username) toggleDrop(exposure);
          break;
        case "0":
          if (username) {
            setOpenPanel("call");
            toggleFormFactor();
          }
          break;
        default: {
          const n = Number(ev.key);
          const c = Number.isInteger(n) && n >= 1 ? candidates[n - 1] : undefined;
          if (c && username) {
            setOpenPanel("call");
            toggleCandidate(c);
          } else handled = false;
        }
      }
      if (handled) ev.preventDefault();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [visible, sampleId, goSample, exposureList, exposure, toggleDrop, candidates, toggleCandidate, username, draft, addCustom, setDraftClamped, openCustom, toggleFormFactor, view, setView, sample, openLabeling]);

  const droppedCount = exposureList.filter((e) => e.status === "rejected").length;
  const keptHint = droppedCount > 0 ? `${exposureList.length - droppedCount} of ${exposureList.length} kept` : undefined;
  const reviewed = samples.filter((s) => sampleStatus(s).kind !== "unreviewed").length;
  const exposureIndex = exposure ? exposureList.findIndex((e) => e.id === exposure.id) + 1 : 0;
  const emptyMessage = !Number.isInteger(experimentId)
    ? "This workspace link is not valid."
    : corpus.isLoading
      ? "Loading samples…"
      : samples.length === 0
        ? ingesting ? "No samples yet. They appear as the ingest reads the data folder." : "No samples in this experiment."
        : "This sample is not in this experiment.";

  const palette = workspacePalette({
    modKey: MOD_KEY,
    sampleName: sample?.name,
    canEdit: Boolean(username && sample),
    nextUnreviewed: () => {
      const nid = nextUnreviewed(visible, sampleId, 1);
      if (nid !== undefined) goSample(nid);
    },
    compare: {
      samples: pickable,
      notes: bufferNotes,
      pick: (rid) => {
        setRefs((r) => addCompare(r, rid));
        if (subRef === undefined) setSubRef(rid);
        setRefsOn(true);
      },
    },
    bridge: {
      samples: bridgePickable,
      pick: (rid) => {
        if (sampleId === undefined) return;
        bridgeSamples.mutate({ sampleId, toSampleId: rid });
        setRefsOn(true);
      },
      bridged: bridged.length > 0,
      leave: leaveBridges,
    },
    candidates: {
      list: allCandidates,
      ticked: assigned,
      tick: (ix) => {
        setOpenPanel("call");
        toggleCandidate(ix);
      },
    },
    labels: {
      open: openLabeling,
      onSample: sampleLabels.map((t) => ({ label: t.label, name: displayLabel(t.label) })),
      recolour: (label) => setColoringLabel(label),
    },
    customIndex: openCustom,
    formFactor: {
      on: isFormFactor,
      toggle: () => {
        setOpenPanel("call");
        toggleFormFactor();
      },
    },
    views: { options: PLOT_VIEWS, current: view, set: setView },
    resetZoom: () => setXDomain(null),
    comparisons: { shown: refsOn, toggle: () => setRefsOn((v) => !v) },
    rescan: ingesting ? undefined : () => triggerScan.mutate(ingestFailed || failedAnalyses > 0),
    problems: problemRows.length > 0 || ingestFailed ? () => setProblemsOpen(true) : undefined,
    removeExperiment:
      ingesting || !username
        ? undefined
        : () => {
            removeExperiment.reset();
            setRemoveOpen(true);
          },
    home: () => navigate("/experiments"),
    theme: { current: theme.theme, choose: theme.choose },
    user: { name: username ?? undefined, switch: clearUsername },
    help: () => setHelpOpen(true),
    samples,
    goSample,
  });

  return (
    <WorkspaceLayout
      theme={theme.theme}
      overlay={
        <>
          {helpOpen && <HelpOverlay sections={HELP} onClose={() => setHelpOpen(false)} />}
          {removeOpen && (
            <RemoveExperimentDialog
              name={experiment.data?.name ?? "this experiment"}
              samples={samples.length}
              removing={removeExperiment.isPending}
              error={removeExperiment.error instanceof Error ? removeExperiment.error.message : null}
              onRemove={() => removeExperiment.mutate(undefined, { onSuccess: () => navigate("/experiments") })}
              onClose={() => setRemoveOpen(false)}
            />
          )}
          {problemsOpen && (
            <ProblemList
              problems={problems.data}
              onPick={(sid) => { setProblemsOpen(false); goSample(sid); }}
              onClose={() => setProblemsOpen(false)}
            />
          )}
          {bigOpen && exposure && (
            <DetectorOverlay
              exposures={exposureList}
              currentId={exposure.id}
              onStep={(dir) => {
                const eid = stepExposure(exposureList, exposure.id, dir);
                if (eid !== undefined) setExposureId(eid);
              }}
              onToggleDrop={username ? () => toggleDrop(exposure) : undefined}
              onClose={() => setBigOpen(false)}
            />
          )}
          {paletteOpen && <Palette groups={palette.groups} pages={palette.pages} onClose={() => setPaletteOpen(false)} />}
        </>
      }
      header={
        <WorkspaceHeader
          experimentName={experiment.data?.name ?? `Experiment ${experimentId}`}
          reviewed={reviewed}
          total={samples.length}
          onHelp={() => setHelpOpen(true)}
        />
      }
      nav={<SampleList samples={visible} currentId={sampleId} filter={filter} onFilter={setFilter} onPick={goSample} bridgeCounts={counts} notice={<ExperimentNotice progress={ingestProgress} ingestFailed={ingestFailed} waiting={waitingForTrace} failed={failedAnalyses} onDetails={() => setProblemsOpen(true)} />} />}
      main={
        sample ? (
          <>
            <SampleHeading
              name={sample.name}
              meta={exposure ? `Frame ${exposureIndex} of ${exposureList.length}` : undefined}
              aside={
                <LabelRow
                  labels={sampleLabels}
                  vocabulary={vocabulary}
                  editable={Boolean(username)}
                  adding={labeling}
                  onOpenAdd={openLabeling}
                  onCloseAdd={() => setLabeling(false)}
                  colors={labelColors.data ?? []}
                  coloring={coloringLabel}
                  onOpenColor={(label) => {
                    setLabeling(false);
                    setColoringLabel(label);
                  }}
                  onCloseColor={() => setColoringLabel(undefined)}
                  onSetColor={(label, color) => setLabelColor.mutate({ label, color })}
                  onAdd={(label) => addLabel.mutate({ sampleId: sample.id, label })}
                  onRemove={(label) => removeLabel.mutate({ sampleId: sample.id, label })}
                />
              }
            />
            <div className="relative px-3.5 py-1">
              <CompareBar
                chips={activeRefs.map((rid) => ({
                  id: rid,
                  name: sampleNames.get(rid) ?? `Sample ${rid}`,
                  color: COMPARE_COLORS[shownRefs.indexOf(rid) % COMPARE_COLORS.length]!,
                  scale: scales.get(rid),
                  bridged: bridged.includes(rid),
                }))}
                subtractId={subRef}
                hidden={!refsOn}
                onSubtract={(rid) => setSubRef((cur) => (cur === rid ? undefined : rid))}
                onRemove={onRemoveChip}
                onAdd={() => setPicking(true)}
              />
              {picking && (
                <ComparePicker
                  samples={pickable}
                  notes={bufferNotes}
                  onPick={(rid) => {
                    setRefs((r) => addCompare(r, rid));
                    if (subRef === undefined) setSubRef(rid);
                    setRefsOn(true);
                  }}
                  onClose={() => setPicking(false)}
                />
              )}
            </div>
            <div className="mx-3.5">
              <PlotCard {...(darkNotice ? { notice: darkNotice } : {})}>
                <WorkspacePlot
                  trace={trace.data}
                  {...(noTraceMessage
                    ? { emptyMessage: noTraceMessage }
                    : trace.isError
                      ? { onRetry: () => void trace.refetch(), retrying: trace.isFetching }
                      : {})}
                  peaks={peaks.data ?? []}
                  combs={combs}
                  compares={compares}
                  xDomain={shownDomain}
                  view={view}
                  onXDomain={setXDomain}
                  {...(draft
                    ? { onAddPeak: snapDraft, onClickPeak: snapDraftToPeak }
                    : canEditPeaks
                      ? { onAddPeak: (q: number) => addPeak.mutate(q), onClickPeak }
                      : {})}
                />
                {difference && xExtent && subRef !== undefined && (
                  <SubtractionStrip
                    q={difference.q}
                    d={difference.d}
                    label={`minus ${sampleNames.get(subRef) ?? `Sample ${subRef}`}`}
                    xDomain={shownDomain}
                    xExtent={xExtent}
                    peakQs={peakQs}
                    view={view}
                  />
                )}
              </PlotCard>
            </div>
            <MainSection title="Exposures" hint={keptHint}>
              <ExposureStrip exposures={exposureList} currentId={exposure?.id} onPick={setExposureId} onToggleDrop={toggleDrop} />
            </MainSection>
          </>
        ) : (
          <div className="px-3.5 py-3">
            <Note>{emptyMessage}</Note>
          </div>
        )
      }
      tools={
        <>
          <Section
            title="Indexing"
            aside={<CustomIndexButton open={draft !== null} disabled={!canCall} onClick={draft ? () => setDraft(null) : openCustom} />}
          >
            <IndexingPlot series={indexSeries} />
          </Section>
          <div className={callOpen ? "flex min-h-24 flex-1 flex-col overflow-y-auto" : "shrink-0"} data-testid="workspace-tools-middle">
            {draft && draftBounds ? (
              <Section title="Custom index">
                <CustomIndexPanel
                  draft={draft}
                  min={draftBounds.min}
                  max={draftBounds.max}
                  landed={draftFit?.landed ?? 0}
                  total={draftFit?.points.length ?? 0}
                  disabled={!canCall}
                  onChange={setDraftClamped}
                  onAdd={addCustom}
                  onCancel={() => setDraft(null)}
                />
              </Section>
            ) : (
              <Section title="Call" collapsible={shortWindow ? { open: callOpen, onOpen: () => setOpenPanel("call") } : undefined}>
                {exposureDropped && <Note>This frame is dropped. Restore it to call a phase.</Note>}
                <CandidateList
                  candidates={candidates}
                  assigned={assigned}
                  disabled={exposureDropped || !username}
                  formFactor={isFormFactor}
                  onToggle={toggleCandidate}
                  onFormFactor={toggleFormFactor}
                  onHover={setHoverId}
                />
              </Section>
            )}
          </div>
          {exposure && (
            <Section
              title="Exposure"
              aside="← →"
              className={shortWindow && exposureOpen ? "min-h-0 flex-1 overflow-y-auto" : exposureOpen ? "min-h-0 shrink overflow-y-auto" : "shrink-0"}
              collapsible={shortWindow ? { open: exposureOpen, onOpen: draft === null ? () => setOpenPanel("exposure") : undefined } : undefined}
            >
              <MetaList rows={exposureRows(exposure)} />
            </Section>
          )}
        </>
      }
    />
  );
}
