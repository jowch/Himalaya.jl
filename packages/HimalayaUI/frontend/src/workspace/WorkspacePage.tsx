import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import type { AssignmentState, CorpusSample, Exposure, IndexEntry, LabelColor, SampleBridge, SampleCallMember, SampleLabel } from "../api";
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
  useSamplesExposures,
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
  useFrameIndices,
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
  useEditOps,
  useCompares,
  usePatchCompares,
  useDetectorSizes,
} from "../queries";
import { detectorGeometry, ringCircles } from "../lib/detectorRings";
import type { DetectorRingsView } from "../ui/DetectorRings";
import { comparesOf, EMPTY_COMPARES } from "../lib/compares";
import { byRecent, recentRank, recordRecent } from "../lib/recents";
import { helpKeyFromEmptyField, inDialog, isEditable, suppressGlobalKeys } from "../lib/keys";
import { useUndoHistory, type UndoEntry, type UndoStep } from "./undo";
import {
  bridgeEntry, sharedAfterBridge, combinedEntry, selectionEntry, settleCallEntry, customIndexEntry, dropEntry, excludeEntry, identify, labelEntry, noPhaseEntry, peakEntry,
  experimentEntry, setupEntry, labelColorEntry, tickEntry, unbridgeEntry, untickEntry, compareEntry, droppedByEdit, type PhaseIdentity,
} from "./undoEntries";
import { showToast } from "../lib/toast";
import { useScanState } from "../hooks/useScanState";
import { useMediaQuery } from "./useMediaQuery";
import { isIngestStage, stageLabel } from "../lib/ingestStages";
import { applyBridge, bridgeCounts, bridgedWith, rankBridgeMates } from "../lib/bridges";
import { bridgeSteps, callStatesSplit, groupKeepChoices, keepChoices, type KeepCall } from "../lib/bridgeCall";
import { bufferRole, buffersFirst, displayLabel, firstUnusedColor, labelColor, labelColorIsDefault, labelsFor, sameLabel, sampleMatchesFilter, labelVocabulary, labelKey } from "../lib/labels";
import { analysisOf, imageFrames, keptExposures, MIN_COMBINED, newestKeptFrame, selectedFrames } from "../lib/exposure";
import { useMemberTraces } from "./useMemberTraces";
import { unsupportedCalls } from "./unsupportedCalls";
import { useCorpusPalette } from "../hooks/useCorpusPalette";
import { allFrames, frameGroups, frameOf } from "./frameGroups";
import { phaseColor } from "../phases";
import { useAppState } from "../state";
import { beamOffWarning, firstUnreviewed, i0BelowZero, imageUrl, nearestKeptExposure, nextUnreviewed, peakClick, rankCandidates, sampleStatus, stepExposure } from "./model";
import { samePeak, type PeakRef } from "../lib/peakRef";
import { clampDomain, compareScale, COMPARES_FULL, MAX_COMPARES, scaledTrace, subtractTrace } from "./compare";
import type { CompareScale } from "./compare";
import { useCompareTraces } from "./useCompareTraces";
import { usePrefetchNeighbours } from "./prefetch";
import { CompareBar } from "../ui/CompareBar";
import { ComparePicker } from "../ui/ComparePicker";
import { SubtractToggle } from "../ui/SubtractToggle";
import { PlotViewToggle } from "../ui/PlotViewToggle";
import { PlotCard } from "../ui/PlotCard";
import { BridgeMenu } from "../ui/BridgeMenu";
import { HelpOverlay, type HelpSection } from "../ui/HelpOverlay";
import { MOD_KEY } from "../lib/modKey";
import { DetectorOverlay } from "../ui/DetectorOverlay";
import { Palette } from "../ui/palette/Palette";
import { keepCallForm, workspacePalette, type PaletteSetup } from "./palette";
import type { PaletteForm } from "../ui/palette/model";
import { useAppTheme } from "../lib/theme";
import { PLOT_VIEWS, stepView, usePlotView } from "./plotView";
import { useSampleListHidden } from "./sampleListHidden";
import { positiveExtent } from "../plot/projection";
import { CustomIndexPanel, NUDGE_STEP, NUDGE_STEP_LARGE } from "../ui/CustomIndexPanel";
import { IndexingPlot, type IndexSeries } from "../ui/IndexingPlot";
import { clampLattice, stepSymmetry, draftCommit, draftPrediction, indexFit, initialDraft, latticeParam, storedFit, type CustomDraft } from "./indexing";
import { latticeBounds, latticeForFirstOrderOnPeak } from "../lib/customIndex";

import { WorkspaceHeader } from "../ui/WorkspaceHeader";
import { Note } from "../ui/Note";
import { WorkspaceLayout } from "../ui/WorkspaceLayout";
import { SampleList } from "../ui/SampleList";
import { ExperimentNotice } from "../ui/ExperimentNotice";
import { ProblemList } from "../ui/ProblemList";
import { RemoveExperimentDialog } from "../ui/RemoveExperimentDialog";
import { COMPARE_COLORS, WorkspacePlot, type CombLine, type CompareLine } from "../ui/WorkspacePlot";
import { ExposureStrip, type ExposureSelection } from "../ui/ExposureStrip";
import { CandidateList } from "../ui/CandidateList";
import { callChoices, cycleChoice, inCallChoices } from "./callTab";
import { callRows, keyedRows } from "./bonnetGroups";
import { MetaList } from "../ui/MetaList";
import { MainSection, SampleHeading, Section } from "../ui/Section";
import { claimedPeaks } from "./claimedPeaks";
import { callCombs } from "./callCombs";
import { exposureKindLabel, exposureTimeLabel } from "./exposureLabels";
import { CustomIndexButton } from "../ui/CustomIndexButton";
import { LabelRow } from "../ui/LabelRow";
import { setupNames } from "../lib/sample/setupNames";
import type { ExperimentPatch } from "../lib/queue/mutators/experiment";
import type { SetupPatch } from "../lib/queue/mutators/setup";

function helpSections(mod: string): HelpSection[] {
  return [
    {
      title: "Move",
      rows: [
        { keys: ["↑", "↓"], what: "Previous or next sample in the list" },
        { keys: ["Z"], what: "Next unreviewed sample" },
        { keys: ["/"], what: "Filter samples by name or label. Press Tab to move through the matches and Enter to open one." },
        { keys: ["["], what: "Hide or show the sample list. While it is hidden, rest the pointer on the left edge to look at it." },
      ],
    },
    {
      title: "Images",
      rows: [
        { keys: ["←", "→"], what: "Previous or next image" },
        { keys: ["X"], what: "Drop the image and show the nearest kept one, or restore a dropped one" },
        { keys: ["Enter"], what: "Select the image to analyze, or take it out of the selection" },
        { keys: ["⇧ click"], what: "Select a range of images to analyze" },
        { keys: ["⌘ click"], what: "Add an image to the selection, or take it out" },
        { keys: ["Space"], what: "Enlarge its detector image; ← →, X and Enter work inside it" },
      ],
    },
    {
      title: "Assignment",
      rows: [
        { keys: ["1–9"], what: "Tick or untick a phase candidate" },
        { keys: ["0"], what: "Mark no phase (form factor), or clear it" },
        { keys: ["Tab", "⇧ Tab"], what: "Step through the phase choices, wrapping at the ends. Press Enter to tick the one in focus." },
        { keys: ["Esc"], what: "Leave the phase choices and let Tab reach every other control until you open another sample" },
        { keys: ["I"], what: "Custom index: press Tab or Shift+Tab to change the symmetry, then drag the lattice, nudge it with ← →, or click a peak to snap to it. Press Enter to add it or Esc to cancel." },
        { keys: ["L"], what: "Label the sample, e.g. Buffer. Press Enter to add it or Esc to close." },
        { keys: [mod, "K"], what: "Open the commands and choose Same sample as… to mark this as the same sample measured in another setup. Space after Tab or an arrow ticks several, and Enter marks them all." },
      ],
    },
    {
      title: "Plot",
      rows: [
        { gesture: "Click", what: "Turn a peak off or on, remove one you added, or add one" },
        { gesture: "Scroll", what: "Zoom q" },
        { keys: [`${mod} 0`], what: "Reset the zoom" },
        { keys: ["V"], what: "Next plot view: log, linear, log-log, q²I. ⇧ V goes back." },
        { keys: ["C"], what: "Compare with another sample" },
        { gesture: `${mod}-click a sample in the list`, what: "Compare with it, or stop comparing" },
        { gesture: "Click a comparison", what: "Hide or show its trace" },
        { gesture: `${mod}-click a comparison`, what: "Plot the current trace minus it, or stop" },
      ],
    },
    {
      title: "Anywhere",
      rows: [
        { keys: [mod, "K"], what: "Commands and go to a sample" },
        { keys: [mod, "Z"], what: "Undo the last edit (calls, drops, peaks, custom index, labels, bridges), even on another sample" },
        { keys: [mod, "⇧", "Z"], what: "Redo it (Ctrl+Y also works)" },
        { keys: ["?"], what: "This list" },
      ],
    },
  ];
}

const HELP = helpSections(MOD_KEY);

const RECENT_LABELS = "label";

const EXPOSURE_LABELS = ["Sample", "File", "Kind", "Acquired", "Exposure", "Flight tube", "Energy", "Detector", "I₀", "I₂"];

function exposureRows(e: Exposure | undefined, owner: string | undefined): [string, string][] {
  if (!e) return EXPOSURE_LABELS.map((label) => [label, "—"]);
  const rows: [string, string][] = [["Sample", owner ?? "—"], ["File", e.filename ?? "—"]];
  rows.push(["Kind", exposureKindLabel(e)]);
  const i0 = e.incident_monitor;
  const i2 = e.transmitted_counts;
  rows.push(["Acquired", e.timestamp ? e.timestamp.replace("T", " ").slice(0, 19) : "—"]);
  rows.push(["Exposure", exposureTimeLabel(e)]);
  rows.push(["Flight tube", e.pipe_length_m != null ? `${Math.round(e.pipe_length_m * 1000)} mm` : "—"]);
  rows.push(["Energy", e.beam_energy_ev != null ? `${(e.beam_energy_ev / 1000).toFixed(3)} keV` : "—"]);
  rows.push(["Detector", e.detector || "—"]);
  rows.push(["I₀", i0 != null ? i0.toLocaleString("en-US", { maximumSignificantDigits: 4 }) : "—"]);
  rows.push(["I₂", i2 != null ? Math.round(i2).toLocaleString("en-US") : "—"]);
  return rows;
}

const SHORT_WINDOW = "(max-height: 860px)";

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
  const unrecognizedSetups = problems.data?.setups?.length ?? 0;
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

  const setups = useMemo(() => experiment.data?.setups ?? [], [experiment.data]);
  const setupOrder = useMemo(() => new Map(setups.map((st, i) => [st.id, i])), [setups]);
  const hiddenSetups = useMemo(() => new Set(setups.filter((st) => !st.included).map((st) => st.id)), [setups]);
  const setupsKnown = experiment.data !== undefined || experiment.isError;
  const inExperiment = useMemo(
    () => (setupsKnown ? (corpus.data ?? []).filter((s) => s.experiment_id === experimentId) : []),
    [setupsKnown, corpus.data, experimentId],
  );
  const samples = useMemo(() => {
    const rank = (s: CorpusSample): number => setupOrder.get(s.setup_id ?? -1) ?? 0;
    return inExperiment
      .filter((s) => s.setup_id == null || !hiddenSetups.has(s.setup_id))
      .map((s, i) => ({ s, i }))
      .sort((a, b) => rank(a.s) - rank(b.s) || a.i - b.i)
      .map(({ s }) => s);
  }, [inExperiment, setupOrder, hiddenSetups]);
  const hiddenSamples = useMemo(
    () => new Set(inExperiment.filter((s) => s.setup_id != null && hiddenSetups.has(s.setup_id)).map((s) => s.id)),
    [inExperiment, hiddenSetups],
  );
  const sampleSetup = useMemo<ReadonlyMap<number, string>>(() => {
    const names = setupNames(setups);
    return new Map(inExperiment.flatMap((s) => {
      const name = s.setup_id == null ? undefined : names?.get(s.setup_id);
      return name ? [[s.id, name] as const] : [];
    }));
  }, [setups, inExperiment]);
  const setupHeaders = useMemo(() => {
    const shown = setups.filter((st) => st.included);
    return shown.length > 1 ? new Map(shown.map((st) => [st.id, st.name])) : undefined;
  }, [setups]);
  const labels = useLabels(experimentId);
  const labelList = useMemo(() => labels.data ?? [], [labels.data]);
  const [recentLabels, setRecentLabels] = useState(() => recentRank(RECENT_LABELS));
  const vocabulary = useMemo(
    () => byRecent(labelVocabulary(labelList), labelKey, recentLabels),
    [labelList, recentLabels],
  );
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

  const bridges = useBridges(experimentId);
  const bridgeList = useMemo(() => bridges.data ?? [], [bridges.data]);
  const shownBridges = useMemo(() => bridgeList.filter((b) => !hiddenSamples.has(b.sample_id)), [bridgeList, hiddenSamples]);
  const bridged = useMemo(() => (sampleId === undefined ? [] : bridgedWith(shownBridges, sampleId)), [shownBridges, sampleId]);
  const exposures = useExposures(sampleId);
  const mateExposures = useSamplesExposures(bridged);
  const [exposureId, setExposureId] = useState<number | undefined>(linkedExposure);
  useEffect(() => setExposureId(linkedExposure), [sampleId, linkedExposure]);
  const ownFrames: Exposure[] = useMemo(() => imageFrames(exposures.data), [exposures.data]);
  const analysis = useMemo(() => analysisOf(exposures.data), [exposures.data]);
  const frameSetupName = useCallback(
    (s: CorpusSample): string => {
      const st = setups.length > 1 && s.setup_id != null ? setups.find((x) => x.id === s.setup_id)?.name : undefined;
      return st ? `${s.name} in ${st}` : s.name;
    },
    [setups],
  );
  const groups = useMemo(() => {
    if (!sample) return [];
    const mates = bridged.flatMap((id, i) => {
      const s = inExperiment.find((x) => x.id === id);
      return s ? [{ sample: s, frames: imageFrames(mateExposures[i]) }] : [];
    });
    return frameGroups({ sample, frames: ownFrames }, mates, (s) => setupOrder.get(s.setup_id ?? -1) ?? 0, frameSetupName);
  }, [sample, bridged, inExperiment, ownFrames, mateExposures, setupOrder, frameSetupName]);
  const exposureList: Exposure[] = useMemo(() => allFrames(groups), [groups]);
  const tile = exposureList.find((e) => e.id === exposureId) ?? analysis.selected[0] ?? newestKeptFrame(ownFrames);
  const exposure = analysis.analyzed;
  const averaging = analysis.selected.length >= MIN_COMBINED;
  const memberIds = useMemo(() => (averaging ? analysis.selected.map((m) => m.id) : []), [averaging, analysis.selected]);
  const members = useMemberTraces(memberIds);
  const tileFrame = frameOf(groups, tile?.id);

  const trace = useTrace(exposure?.id);
  const noTrace = trace.error instanceof ApiError && trace.error.status === 404;
  const i0Hint = i0BelowZero(exposure);
  const darkNotice = exposure?.beam_off && trace.data
    ? `Beam likely off: the detector recorded no counts, so the beam was off or the shutter closed.${i0Hint ? " I₀ also reads below zero." : ""}`
    : undefined;
  const tracePending = problemRows.some((p) => p.problem === "no_trace" && p.exposure_id === exposure?.id);
  const selectionMessage = !exposures.data || ownFrames.length === 0
    ? undefined
    : analysis.selected.length === 0
      ? "Select an image to analyze."
      : exposure === undefined
        ? analysis.selected.some((e) => e.trace_hash === null)
          ? "Waiting for the selected images' traces."
          : "Averaging the selected images…"
        : undefined;
  const noTraceMessage = selectionMessage ?? (!noTrace
    ? undefined
    : tracePending
      ? "Waiting for this image's trace."
      : averaging
        ? "Couldn't average the selected images. They may not share a q range."
        : `No integrated trace for this image.${i0Hint ? " I₀ reads below zero, so the beam was likely off." : ""}`);
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
  const candidates = useMemo(() => rankCandidates(indices.data), [indices.data]);
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
  const offeredPhases = useMemo(() => new Set(candidates.map((c) => c.phase)), [candidates]);
  const unsupportedFrames = useMemo(
    () => [...new Set(otherMembers.filter((m) => !offeredPhases.has(m.phase)).map((m) => m.exposure_id))],
    [otherMembers, offeredPhases],
  );
  const unsupportedLattices = useFrameIndices(unsupportedFrames);
  const unsupported = useMemo(
    () => unsupportedCalls(otherMembers, candidates, unsupportedLattices),
    [otherMembers, candidates, unsupportedLattices],
  );
  const unsupportedIds = useMemo(() => new Set(unsupported.map((c) => c.id)), [unsupported]);
  const choices = useMemo(() => [...candidates, ...unsupported], [candidates, unsupported]);
  const assigned = useMemo(() => {
    const here = new Set(assignment.data?.members ?? []);
    const phases = new Set(otherMembers.map((m) => m.phase));
    for (const c of choices) if (phases.has(c.phase)) here.add(c.id);
    return here;
  }, [assignment.data, otherMembers, choices]);
  const rows = useMemo(
    () => [...callRows(candidates), ...unsupported.map((candidate) => ({ kind: "one" as const, candidate }))],
    [candidates, unsupported],
  );
  const keyed = useMemo(() => keyedRows(rows), [rows]);
  const [hoverIds, setHoverIds] = useState<readonly number[] | undefined>(undefined);
  const hoverCandidates = useCallback((ids: readonly number[] | undefined) => {
    setHoverIds((prev) => (prev?.length === ids?.length && (prev ?? []).every((id, i) => id === ids?.[i]) ? prev : ids));
  }, []);

  const addPhase = useAddAssignmentPhase(exposure?.id ?? 0);
  const removePhase = useRemoveAssignmentPhase(exposure?.id ?? 0);
  const removeFramePhase = useRemoveFramePhase();
  const setCallState = useSetSampleCallState(sampleId ?? 0);
  const isFormFactor = assignment.data?.state === "form_factor";
  const setStatus = useSetExposureStatus(sampleId ?? 0);

  const canEditPeaks = Boolean(username && exposure);
  const addPeak = useAddPeak(exposure?.id ?? 0);
  const removePeak = useRemovePeak(exposure?.id ?? 0);
  const setExcluded = useSetPeakExcluded(exposure?.id ?? 0);
  const ops = useEditOps(experimentId);
  const history = useUndoHistory(`${experimentId}:${username ?? ""}`, ops.whenSaved);
  const frameName = useCallback(
    (eid: number) => {
      if (averaging && eid === exposure?.id) return `the average of ${analysis.selected.length} images of ${sample?.name ?? "this sample"}`;
      const at = frameOf(groups, eid);
      const owner = at && at.group.sampleId !== sample?.id ? inExperiment.find((s) => s.id === at.group.sampleId)?.name : sample?.name;
      return `image ${(at?.index ?? -1) + 1} of ${owner ?? "this sample"}`;
    },
    [groups, sample?.name, sample?.id, inExperiment, averaging, exposure?.id, analysis.selected.length],
  );

  const [xDomain, setXDomain] = useState<[number, number] | null>(null);
  const [view, setView] = usePlotView();
  const [sampleListHidden, setSampleListHidden] = useSampleListHidden();
  const comparesQuery = useCompares(experimentId);
  const patchCompares = usePatchCompares(experimentId);
  const saved = useMemo(() => comparesOf(comparesQuery.data ?? EMPTY_COMPARES, sampleId ?? -1), [comparesQuery.data, sampleId]);
  const refs = saved.refs;
  const subChoice = saved.subtraction;
  const [refsOn, setRefsOn] = useState(true);
  const [shownOverride, setShownOverride] = useState<ReadonlyMap<string, boolean>>(() => new Map());
  const [picking, setPicking] = useState(false);
  const [helpOpen, setHelpOpen] = useState(false);
  const [bigOpen, setBigOpen] = useState(false);
  const [paletteOpen, setPaletteOpen] = useState(false);
  const paletteCorpus = useCorpusPalette();
  const elsewhere = useMemo(
    () => ({
      experiments: paletteCorpus.experiments.filter((e) => e.id !== experimentId),
      samples: paletteCorpus.samples.filter((s) => s.experimentId !== experimentId),
    }),
    [paletteCorpus.experiments, paletteCorpus.samples, experimentId],
  );
  const theme = useAppTheme();
  const [draft, setDraft] = useState<CustomDraft | null>(null);
  const customPanel = useRef<HTMLDivElement>(null);
  useEffect(() => setDraft(null), [exposure?.id]);
  useEffect(() => setHoverIds(undefined), [exposure?.id]);
  const commitCustom = useCommitCustomIndex(exposure?.id ?? 0);
  const bridgedNames = useMemo(() => {
    const names = new Map(inExperiment.map((s) => [s.id, s.name]));
    return bridged.map((r) => names.get(r) ?? `Sample ${r}`);
  }, [bridged, inExperiment]);
  const bridgeMates = useMemo(
    () => bridged.map((r) => {
      const s = inExperiment.find((x) => x.id === r);
      return { id: r, title: s ? frameSetupName(s) : `Sample ${r}` };
    }),
    [bridged, inExperiment, frameSetupName],
  );
  const counts = useMemo(() => bridgeCounts(shownBridges), [shownBridges]);
  const bridgeSamples = useBridgeSamples(experimentId);
  const unbridgeSample = useUnbridgeSample(experimentId);
  const shownRefs = useMemo(() => refs.filter((r) => !hiddenSamples.has(r)), [refs, hiddenSamples]);
  const subRef = subChoice !== undefined && shownRefs.includes(subChoice) ? subChoice : undefined;
  const hiddenRefs = useMemo<ReadonlySet<number>>(
    () =>
      new Set(
        shownRefs.filter((r) => {
          const shown = shownOverride.get(`${sampleId}:${r}`);
          if (shown !== undefined) return !shown;
          const by = saved.addedBy.get(r);
          return by !== undefined && by !== (username ?? null);
        }),
      ),
    [shownRefs, shownOverride, sampleId, saved.addedBy, username],
  );
  const setShown = useCallback(
    (rid: number, shown: boolean) => setShownOverride((m) => new Map(m).set(`${sampleId}:${rid}`, shown)),
    [sampleId],
  );
  const fetchRefs = useMemo(() => shownRefs.filter((r) => !hiddenRefs.has(r) || r === subRef), [shownRefs, hiddenRefs, subRef]);
  const compareTraces = useCompareTraces(fetchRefs);

  const callOnFrame = useCallback(
    (eid: number): PhaseIdentity[] => {
      if (assignment.data?.exposure_id !== eid) return [];
      const members = new Set(assignment.data.members);
      return (indices.data ?? [])
        .filter((ix) => members.has(ix.id))
        .map((ix) => ({ exposureId: eid, phase: ix.phase, basis: ix.basis, kind: ix.kind }));
    },
    [assignment.data, indices.data],
  );

  const onClickPeak = useCallback(
    (ref: PeakRef) => {
      const action = peakClick(peaks.data ?? [], ref);
      if (!action || !exposure || sampleId === undefined) return;
      const eid = exposure.id;
      const q = (peaks.data ?? []).find((pk) => samePeak(pk, ref))?.q;
      const before = callOnFrame(eid);
      if (action.kind === "remove") {
        if (ref.id < 0) return;
        removePeak.mutate(ref.id);
        if (q !== undefined) {
          history.record(peakEntry(ops, {
            label: `removed a peak at q = ${q.toFixed(4)} on ${frameName(eid)}`, sampleId, exposureId: eid, q, added: false,
            dropped: droppedByEdit(ops, before),
          }));
        }
        return;
      }
      setExcluded.mutate({ peakId: ref.id, excluded: action.excluded });
      if (q !== undefined) {
        history.record(excludeEntry(ops, {
          label: `${action.excluded ? "turned off" : "turned on"} a peak on ${frameName(eid)}`, sampleId, exposureId: eid, q,
          excluded: action.excluded, dropped: droppedByEdit(ops, before),
        }));
      }
    },
    [peaks.data, removePeak, setExcluded, exposure, sampleId, history, ops, frameName, callOnFrame],
  );

  const addPeakAt = useCallback(
    (q: number) => {
      if (!exposure || sampleId === undefined) return;
      const eid = exposure.id;
      addPeak.mutate(q);
      history.record(peakEntry(ops, {
        label: `added a peak at q = ${q.toFixed(4)} on ${frameName(eid)}`, sampleId, exposureId: eid, q, added: true,
        dropped: droppedByEdit(ops, callOnFrame(eid)),
      }));
    },
    [exposure, sampleId, addPeak, history, ops, frameName, callOnFrame],
  );


  const callEdit = useCallback(
    (ix: IndexEntry, wasFormFactor: boolean): UndoEntry | undefined => {
      if (!exposure || sampleId === undefined) return undefined;
      const eid = exposure.id;
      const where = sample?.name ?? "this sample";
      if (assignment.data?.members.includes(ix.id)) {
        removePhase.mutate(ix.id);
        return untickEntry(ops, {
          label: `unticked ${ix.phase} on ${where}`, sampleId, exposureId: eid,
          ids: Promise.resolve([{ exposureId: eid, phase: ix.phase, basis: ix.basis, kind: ix.kind }]),
        });
      }
      const elsewhere = otherMembers.filter((m) => m.phase === ix.phase);
      if (elsewhere.length > 0) {
        const ids = identify(ops, elsewhere.map((m) => ({ exposureId: m.exposure_id, indexId: m.index_id })));
        for (const m of elsewhere) removeFramePhase.mutate(m.exposure_id, m.index_id);
        return untickEntry(ops, { label: `unticked ${ix.phase} on ${where}`, sampleId, exposureId: eid, ids });
      }
      if (exposure.status === "rejected") return undefined;
      addPhase.mutate(ix.id);
      return tickEntry(ops, {
        label: `ticked ${ix.phase} on ${where}`, sampleId, exposureId: eid,
        id: { exposureId: eid, phase: ix.phase, basis: ix.basis, kind: ix.kind }, wasFormFactor,
      });
    },
    [exposure, sampleId, sample?.name, assignment.data, otherMembers, addPhase, removePhase, removeFramePhase, ops],
  );

  const toggleCandidate = useCallback(
    (ix: IndexEntry) => {
      const entry = callEdit(ix, isFormFactor);
      if (entry) history.record(entry);
    },
    [callEdit, isFormFactor, history],
  );

  const toggleGroup = useCallback(
    (members: readonly IndexEntry[]) => {
      if (!exposure || sampleId === undefined) return;
      const allIn = members.every((m) => assigned.has(m.id));
      const parts = (allIn ? members : members.filter((m) => !assigned.has(m.id)))
        .map((m, i) => callEdit(m, isFormFactor && i === 0))
        .filter((e): e is UndoEntry => e !== undefined);
      if (parts.length === 0) return;
      const phases = members.map((m) => m.phase).join(" + ");
      history.record(combinedEntry({
        label: `${allIn ? "unticked" : "ticked"} ${phases} on ${sample?.name ?? "this sample"}`, sampleId, exposureId: exposure.id,
      }, parts));
    },
    [exposure, sampleId, sample?.name, assigned, callEdit, isFormFactor, history],
  );

  const toggleFormFactor = useCallback(() => {
    if (!exposure || exposure.status === "rejected" || sampleId === undefined) return;
    const on = !isFormFactor;
    const eid = exposure.id;
    const ticked = on
      ? identify(ops, [
          ...(assignment.data?.members ?? []).map((indexId) => ({ exposureId: eid, indexId })),
          ...otherMembers.map((m) => ({ exposureId: m.exposure_id, indexId: m.index_id })),
        ])
      : Promise.resolve([]);
    setCallState.mutate(on ? "form_factor" : "null");
    history.record(noPhaseEntry(ops, {
      label: `${on ? "marked" : "cleared"} no phase on ${sample?.name ?? "this sample"}`, sampleId, exposureId: eid, on, ticked,
    }));
  }, [exposure, sampleId, sample?.name, isFormFactor, assignment.data, otherMembers, setCallState, history, ops]);

  const toggleDrop = useCallback(
    (e: Exposure) => {
      const prev = e.status;
      const next = prev === "rejected" ? null : "rejected";
      if (next === null) setExposureId(e.id);
      else if (e.id === tile?.id) setExposureId(nearestKeptExposure(exposureList, e.id) ?? e.id);
      if (sampleId === undefined) return;
      const owner = e.sample_id;
      if (owner === sampleId) setStatus.mutate({ exposureId: e.id, status: next });
      else ops.setStatus(owner, e.id, next);
      history.record(dropEntry(ops, {
        label: `${next ? "dropped" : "restored"} ${frameName(e.id)}`, sampleId, exposureId: e.id, owner, prev, next,
        ...(next === "rejected" ? { selected: e.selected !== false } : {}),
      }));
    },
    [setStatus, sampleId, history, ops, frameName, exposureList, tile],
  );
  const anchors = useRef(new Map<number, number>());
  useEffect(() => { anchors.current = new Map(); }, [sampleId]);
  const groupSelections = useMemo(() => {
    if (sampleId === undefined) return new Map<number, ExposureSelection>();
    return new Map(groups.map((g) => {
      const owner = g.sampleId;
      const ownerName = (owner === sampleId ? sample?.name : inExperiment.find((x) => x.id === owner)?.name) ?? "this sample";
      const kept = keptExposures(g.exposures).map((e) => e.id);
      const prev = selectedFrames(g.exposures).map((e) => e.id);
      const selectedIds = new Set(prev);
      const apply = (next: number[], label: string): void => {
        ops.setSelection(owner, next);
        history.record(selectionEntry(ops, { label, sampleId, owner, exposureId: undefined, prev, next }));
      };
      const selection: ExposureSelection = {
        selectedIds,
        selectableIds: new Set(kept),
        onToggle: (eid) => {
          anchors.current.set(owner, eid);
          const on = selectedIds.has(eid);
          apply(on ? prev.filter((id) => id !== eid) : [...prev, eid], `${on ? "deselected" : "selected"} ${frameName(eid)}`);
        },
        onRange: (eid) => {
          const to = kept.indexOf(eid);
          if (to < 0) return;
          const at = kept.indexOf(anchors.current.get(owner) ?? prev[0] ?? eid);
          const from = at < 0 ? to : at;
          const span = kept.slice(Math.min(from, to), Math.max(from, to) + 1);
          anchors.current.set(owner, eid);
          apply([...new Set([...prev, ...span])], `selected ${span.length} images of ${ownerName}`);
        },
        onAll: () => apply(kept, `selected every image of ${ownerName}`),
        onNone: () => apply([], `cleared the selection of ${ownerName}`),
      };
      return [owner, selection] as const;
    }));
  }, [groups, sampleId, sample?.name, inExperiment, ops, history, frameName]);

  const tileSelection = tile && username ? groupSelections.get(tile.sample_id) : undefined;
  const tileSelectable = tile !== undefined && (tileSelection?.selectableIds.has(tile.id) ?? false);
  const tileSelected = tile !== undefined && (tileSelection?.selectedIds.has(tile.id) ?? false);
  const toggleTileSelection = useCallback(() => {
    if (tile && tileSelectable) tileSelection?.onToggle(tile.id);
  }, [tile, tileSelectable, tileSelection]);

  const draftQs = useMemo(() => (draft ? draftPrediction(draft) : []), [draft]);
  const combs: CombLine[] = useMemo(
    () => callCombs(choices, hoverIds, assigned, draft ? { phase: draft.phase, qs: draftQs } : null),
    [choices, hoverIds, assigned, draft, draftQs],
  );

  const ringSets = useMemo(() => [
    ...choices.filter((c) => assigned.has(c.id)).map((c) => ({ color: phaseColor(c.phase), qs: c.predicted_q })),
    ...(draft ? [{ color: phaseColor(draft.phase), qs: draftQs }] : []),
  ], [draft, draftQs, choices, assigned]);
  const setupOfSample = useMemo(
    () => new Map(inExperiment.map((s) => [s.id, setups.find((st) => st.id === s.setup_id)] as const)),
    [inExperiment, setups],
  );
  const qUnits = setupOfSample.get(sampleId ?? -1)?.q_units;
  const geometryOf = useCallback((e: Exposure) => detectorGeometry(setupOfSample.get(e.sample_id ?? -1)), [setupOfSample]);
  const sizeUrls = useMemo(
    () => exposureList.flatMap((e) => {
      const url = geometryOf(e) ? imageUrl(e, true) : null;
      return url ? [url] : [];
    }),
    [exposureList, geometryOf],
  );
  const detectorSizes = useDetectorSizes(sizeUrls);
  const ringsOf = useCallback(
    (e: Exposure): DetectorRingsView | null => {
      if (ringSets.length === 0) return null;
      const g = geometryOf(e);
      const url = imageUrl(e, true);
      const size = url ? detectorSizes.get(url) : undefined;
      const centre = g ? ringCircles([], qUnits, g, size) : null;
      if (!g || !size || !centre) return null;
      const rings = ringSets.flatMap((set) => (ringCircles(set.qs, qUnits, g, size)?.radii ?? []).map((r) => ({ r, color: set.color })));
      return { size, cx: centre.cx, cy: centre.cy, side: centre.side, rings };
    },
    [ringSets, geometryOf, detectorSizes, qUnits],
  );

  const highlight = useMemo(
    () => (draft ? null : claimedPeaks(candidates.filter((c) => hoverIds?.includes(c.id)), peaks.data ?? [])),
    [draft, candidates, hoverIds, peaks.data],
  );

  const sampleNames = useMemo(() => new Map(samples.map((s) => [s.id, s.name])), [samples]);
  const nameOf = useCallback((rid: number) => sampleNames.get(rid) ?? `sample ${rid}`, [sampleNames]);
  const removeCompare = useCallback(
    (rid: number) => {
      if (sampleId === undefined) return;
      patchCompares.mutate({ sampleId, change: { remove: [rid] } });
      history.record(compareEntry(ops, {
        label: `stopped comparing ${sample?.name ?? "this sample"} with ${nameOf(rid)}`, sampleId, exposureId: undefined,
        edit: { kind: "remove", ref: rid, subtracted: subChoice === rid },
      }));
    },
    [sampleId, patchCompares, history, ops, sample?.name, nameOf, subChoice],
  );
  const activeRefs = useMemo(() => shownRefs.filter((r) => r !== sampleId), [shownRefs, sampleId]);
  const comparesFull = activeRefs.length >= MAX_COMPARES;
  const compareColors = useMemo(
    () => new Map(activeRefs.map((rid, i) => [rid, COMPARE_COLORS[i % COMPARE_COLORS.length]!] as const)),
    [activeRefs],
  );
  const mainTrace = trace.data;
  const peakQs = useMemo(() => (peaks.data ?? []).map((p) => p.q), [peaks.data]);
  const scales = useMemo(() => {
    const out = new Map<number, CompareScale>();
    if (!mainTrace) return out;
    const setupOf = (sid: number | undefined) => inExperiment.find((s) => s.id === sid)?.setup_id ?? null;
    const own = setupOf(sampleId);
    const sameKnownSetup = (rid: number) => own !== null && setupOf(rid) === own;
    compareTraces.forEach((t, rid) => out.set(rid, compareScale(mainTrace, t, peakQs, sameKnownSetup(rid))));
    return out;
  }, [mainTrace, compareTraces, peakQs, inExperiment, sampleId]);
  const compares: CompareLine[] = useMemo(() => {
    if (!mainTrace || !refsOn) return [];
    return activeRefs.flatMap((rid, i) => {
      const t = compareTraces.get(rid);
      const scale = scales.get(rid)?.scale;
      if (rid === sampleId || hiddenRefs.has(rid) || !t || scale === undefined) return [];
      const sc = scaledTrace(t, scale);
      return [{ id: rid, name: sampleNames.get(rid) ?? `Sample ${rid}`, q: sc.q, I: sc.I, color: COMPARE_COLORS[i % COMPARE_COLORS.length]!, emphasis: rid === subRef }];
    });
  }, [mainTrace, refsOn, activeRefs, hiddenRefs, compareTraces, scales, sampleId, sampleNames, subRef]);
  const subTrace = subRef !== undefined && subRef !== sampleId && refsOn ? compareTraces.get(subRef) : undefined;
  const subScale = subRef !== undefined ? scales.get(subRef)?.scale : undefined;
  const difference = useMemo(
    () => (mainTrace && subTrace && subScale !== undefined ? subtractTrace(mainTrace, subTrace, subScale) : null),
    [mainTrace, subTrace, subScale],
  );
  const [minusSamples, setMinusSamples] = useState<ReadonlySet<number>>(() => new Set());
  const setMinus = useCallback(
    (on: boolean) => {
      if (sampleId === undefined) return;
      setMinusSamples((prev) => {
        const next = new Set(prev);
        if (on) next.add(sampleId);
        else next.delete(sampleId);
        return next;
      });
    },
    [sampleId],
  );
  const minusOn = difference !== null && subRef !== undefined && sampleId !== undefined && minusSamples.has(sampleId);
  const subName = subRef !== undefined ? (sampleNames.get(subRef) ?? `Sample ${subRef}`) : "";
  const shownCompareSet = useMemo<ReadonlySet<number>>(
    () => (!refsOn ? new Set(shownRefs) : minusOn ? new Set(shownRefs.filter((r) => r !== subRef)) : hiddenRefs),
    [refsOn, minusOn, shownRefs, subRef, hiddenRefs],
  );
  const xExtent = useMemo<[number, number] | null>(
    () => (mainTrace && mainTrace.q.length > 1 ? positiveExtent(mainTrace.q) : null),
    [mainTrace],
  );
  const shownDomain = clampDomain(xDomain, xExtent);
  const recentKind = `comparison:${experimentId}`;
  const [recentCompares, setRecentCompares] = useState<string[]>(() => recentRank(recentKind));
  const pickable = useMemo(
    () => byRecent(buffersFirst(samples.filter((s) => s.id !== sampleId && !shownRefs.includes(s.id)), labelList), (s) => String(s.id), recentCompares),
    [samples, sampleId, shownRefs, labelList, recentCompares],
  );
  const bufferNotes = useMemo(() => {
    const notes = new Map<number, string>();
    for (const s of samples) {
      const role = bufferRole(labelList, s.id);
      if (role) notes.set(s.id, role === "confirmed" ? "Buffer" : "Suggested buffer");
    }
    return notes;
  }, [samples, labelList]);
  const bridgePickable = useMemo(() => {
    const others = samples.filter((s) => s.id !== sampleId && !bridged.includes(s.id));
    return sample ? rankBridgeMates(sample, others) : others;
  }, [samples, sample, sampleId, bridged]);
  useEffect(() => {
    setLabeling(false);
    setColoringLabel(undefined);
  }, [sampleId]);
  const sampleLabels = useMemo(() => (sampleId === undefined ? [] : labelsFor(labelList, sampleId)), [labelList, sampleId]);
  const paletteLabels = useMemo(
    () =>
      vocabulary.map((name) => {
        const mine = sampleLabels.find((l) => sameLabel(l.label, name));
        return {
          label: name,
          name: displayLabel(name),
          color: labelColor(labelColors.data ?? [], name),
          uses: labelList.filter((l) => l.state === "confirmed" && sameLabel(l.label, name)).length,
          onSample: mine?.state === "confirmed" ? ("confirmed" as const) : mine ? ("suggested" as const) : undefined,
        };
      }),
    [vocabulary, sampleLabels, labelColors.data, labelList],
  );
  const openLabeling = useCallback(() => {
    if (!username || sampleId === undefined) return;
    setColoringLabel(undefined);
    setLabeling(true);
  }, [username, sampleId]);
  const unbridgeOne = useCallback(
    (target: number) => {
      if (!username || sampleId === undefined) return;
      unbridgeSample.mutate({ sampleId: target });
      history.record(unbridgeEntry(ops, {
        label: `marked ${sampleNames.get(target) ?? "a sample"} as not the same sample`, sampleId, exposureId: undefined, target, bridges: bridgeList,
      }));
    },
    [username, sampleId, unbridgeSample, history, ops, sampleNames, bridgeList],
  );
  const leaveBridges = useCallback(() => {
    if (sampleId !== undefined && bridged.length > 0) unbridgeOne(sampleId);
  }, [sampleId, bridged, unbridgeOne]);
  const bridgeTo = useCallback(
    (rid: number, keep?: KeepCall) => {
      if (sampleId === undefined) return;
      const states = {
        from: sample?.assignment_state ?? "indexed",
        to: samples.find((s) => s.id === rid)?.assignment_state ?? "indexed",
      };
      const send = (lost: readonly SampleCallMember[]): void => {
        bridgeSamples.mutate(keep === undefined ? { sampleId, toSampleId: rid } : { sampleId, toSampleId: rid, keepCall: keep });
        history.record(bridgeEntry(ops, {
          label: `marked ${sample?.name ?? "this sample"} as the same sample as ${sampleNames.get(rid) ?? "a sample"}`, sampleId, exposureId: undefined,
          from: sampleId, to: rid, bridges: bridgeList, labels: labelList, keep, states, lost,
        }));
      };
      if (typeof keep !== "number") return send([]);
      ops.read.call(keep === sampleId ? rid : sampleId).then((c) => send(c.members), () => send([]));
    },
    [sampleId, sample?.name, sample?.assignment_state, samples, bridgeSamples, history, ops, sampleNames, bridgeList, labelList],
  );
  const bridgeMany = useCallback(
    async (ids: readonly number[], keep?: KeepCall) => {
      if (sampleId === undefined || !sample) return;
      const mates = ids.flatMap((id) => samples.filter((s) => s.id === id));
      let bridges: SampleBridge[] = [...bridgeList];
      let labels: SampleLabel[] = [...labelList];
      let groupState: AssignmentState = sample.assignment_state ?? "indexed";
      const parts = [];
      for (const step of bridgeSteps(sample, mates, keep)) {
        if (bridgedWith(bridges, sampleId).includes(step.to)) continue;
        const mate = mates.find((m) => m.id === step.to)!;
        const states = { from: groupState, to: mate.assignment_state ?? "indexed" };
        if (typeof step.keep === "number" && parts.length > 0) await ops.whenSaved();
        const lost = typeof step.keep === "number"
          ? await ops.read.call(step.keep === sampleId ? step.to : sampleId).then((c) => c.members, () => [])
          : [];
        bridgeSamples.mutate(step.keep === undefined ? { sampleId, toSampleId: step.to } : { sampleId, toSampleId: step.to, keepCall: step.keep });
        parts.push(bridgeEntry(ops, { label: "", sampleId, exposureId: undefined, from: sampleId, to: step.to, bridges, labels, keep: step.keep, states, lost }));
        bridges = applyBridge(bridges, sampleId, step.to);
        labels = sharedAfterBridge(labels, [sampleId, ...bridgedWith(bridges, sampleId)]);
        if (step.keep === step.to || (groupState === "null" && states.to !== "null")) groupState = states.to;
      }
      history.record(combinedEntry({
        label: `marked ${sample.name} as the same sample as ${mates.map((m) => sampleNames.get(m.id) ?? m.name).join(", ")}`, sampleId, exposureId: undefined,
      }, parts));
    },
    [sampleId, sample, samples, bridgeList, labelList, ops, bridgeSamples, history, sampleNames],
  );
  const nameWithSetup = useCallback(
    (s: Pick<CorpusSample, "id" | "name">): string => {
      const setup = sampleSetup.get(s.id);
      return setup ? `${s.name} in ${setup}` : s.name;
    },
    [sampleSetup],
  );
  const bridgeChoices = useCallback(
    (rid: number) => {
      const there = samples.find((s) => s.id === rid);
      return sample && there ? keepChoices(sample, there, nameWithSetup) : undefined;
    },
    [sample, samples, nameWithSetup],
  );
  const unsettledMate = useMemo(() => {
    if (!sample) return undefined;
    return bridged.map((id) => samples.find((s) => s.id === id)).find((m) => m !== undefined && callStatesSplit(sample, m));
  }, [sample, samples, bridged]);
  const [paletteStart, setPaletteStart] = useState<PaletteForm | undefined>(undefined);
  const settlePrompted = useRef(new Set<number>());
  useEffect(() => {
    if (!username || !sample || !unsettledMate || settlePrompted.current.has(sample.id)) return;
    const choices = keepChoices(sample, unsettledMate, nameWithSetup);
    if (!choices) return;
    settlePrompted.current.add(sample.id);
    const here = sample.id;
    const mate = unsettledMate.id;
    const stateOf = new Map([[here, sample.assignment_state ?? "indexed"], [mate, unsettledMate.assignment_state ?? "indexed"]]);
    const settle = (keep: KeepCall): void => {
      if (typeof keep !== "number") return;
      const state = stateOf.get(keep) ?? "indexed";
      const record = (lost: readonly SampleCallMember[]): void => {
        ops.bridge(here, mate, keep);
        history.record(settleCallEntry(ops, {
          label: `kept ${nameWithSetup(keep === here ? sample : unsettledMate)}'s call for the same sample`, sampleId: here, exposureId: undefined,
          from: here, to: mate, keep, lost,
        }));
      };
      if (state !== "form_factor") return record([]);
      ops.read.call(here).then((c) => record(c.members), () => record([]));
    };
    setPaletteStart(keepCallForm(`Calls differ: same sample as ${nameWithSetup(unsettledMate)}`, choices, settle, "keep"));
    setPaletteOpen(true);
  }, [username, sample, unsettledMate, nameWithSetup, ops, history]);
  const labelEdit = useCallback(
    (sid: number, tag: string, added: boolean) => {
      const prev = labelList.find((l) => l.sample_id === sid && sameLabel(l.label, tag))?.state;
      (added ? addLabel : removeLabel).mutate({ sampleId: sid, label: tag });
      if (added) {
        recordRecent(RECENT_LABELS, labelKey(tag));
        setRecentLabels(recentRank(RECENT_LABELS));
      }
      const verb = prev === "suggested" ? (added ? "confirmed the suggested" : "dismissed the suggested") : added ? "added the" : "removed the";
      return labelEntry(ops, {
        label: `${verb} label ${displayLabel(tag)} on ${sample?.name ?? "this sample"}`,
        sampleId: sid, exposureId: undefined, tag, added, prev,
      });
    },
    [sample?.name, labelList, addLabel, removeLabel, ops],
  );
  const labelRecolor = useCallback(
    (sid: number, tag: string, color: LabelColor | null) => {
      const prev = (labelColors.data ?? []).find((t) => sameLabel(t.label, tag))?.color ?? null;
      setLabelColor.mutate({ label: tag, color });
      return labelColorEntry(ops, {
        label: `recoloured the label ${displayLabel(tag)}`, sampleId: sid, exposureId: undefined, tag, prev, next: color,
      });
    },
    [labelColors.data, setLabelColor, ops],
  );
  const editLabel = useCallback(
    (tag: string, added: boolean) => {
      if (sampleId !== undefined) history.record(labelEdit(sampleId, tag, added));
    },
    [sampleId, history, labelEdit],
  );
  const recolorLabel = useCallback(
    (tag: string, color: LabelColor | null) => {
      if (sampleId !== undefined) history.record(labelRecolor(sampleId, tag, color));
    },
    [sampleId, history, labelRecolor],
  );
  const editExperiment = useCallback(
    (next: ExperimentPatch) => {
      const cur = experiment.data;
      if (!cur) return;
      const prev: ExperimentPatch = {};
      if (next.name !== undefined) prev.name = cur.name ?? "";
      if (next.description !== undefined) prev.description = cur.description ?? null;
      ops.patchExperiment(next);
      if (sampleId === undefined) return;
      const label = next.name !== undefined ? `renamed the experiment to ${next.name}` : "changed the experiment description";
      history.record(experimentEntry(ops, { label, sampleId, exposureId: undefined, prev, next }));
    },
    [experiment.data, ops, sampleId, history],
  );
  const paletteSetups: PaletteSetup[] = useMemo(
    () => setups.map((st) => ({ id: st.id, folder: st.folder, name: st.name, included: st.included })),
    [setups],
  );
  const editSetup = useCallback(
    (setupId: number, next: SetupPatch) => {
      const cur = setups.find((st) => st.id === setupId);
      if (!cur) return;
      const prev: SetupPatch = {};
      if (next.name !== undefined) prev.name = cur.display_name;
      if (next.included !== undefined) prev.included = cur.included;
      ops.patchSetup(setupId, next);
      if (sampleId === undefined) return;
      const label = next.included === false
        ? `hid setup ${cur.name}`
        : next.included === true
          ? `showed setup ${cur.name}`
          : `renamed setup ${cur.name} to ${next.name ?? cur.folder}`;
      history.record(setupEntry(ops, { label, sampleId, exposureId: undefined, setupId, prev, next }));
    },
    [setups, ops, sampleId, history],
  );
  const createLabel = useCallback(
    (name: string, color: LabelColor) => {
      if (sampleId === undefined) return;
      const known = labelColors.data ?? [];
      const recolor = labelColorIsDefault(known, name) || labelColor(known, name) !== color;
      const parts = [...(recolor ? [labelRecolor(sampleId, name, color)] : []), labelEdit(sampleId, name, true)];
      history.record(combinedEntry({
        label: `created the label ${displayLabel(name)} on ${sample?.name ?? "this sample"}`, sampleId, exposureId: undefined,
      }, parts));
    },
    [sampleId, sample?.name, labelColors.data, labelRecolor, labelEdit, history],
  );

  const activePeakQs = useMemo(() => (peaks.data ?? []).filter((p) => !p.excluded).map((p) => p.q), [peaks.data]);
  const canCall = Boolean(username && exposure);
  const draftBounds = draft ? latticeBounds(draft.phase, xExtent) : null;
  const shortWindow = useMediaQuery(SHORT_WINDOW);
  const [openPanel, setOpenPanel] = useState<"call" | "exposure">("call");
  const callOpen = !shortWindow || openPanel === "call" || draft !== null;
  const exposureOpen = !shortWindow || (openPanel === "exposure" && draft === null);
  const tabFree = useRef(false);
  const [pendingTab, setPendingTab] = useState<"first" | "last" | null>(null);
  useEffect(() => {
    tabFree.current = false;
  }, [sampleId]);
  useEffect(() => {
    if (!pendingTab || !callOpen) return;
    const choices = callChoices();
    (pendingTab === "last" ? choices[choices.length - 1] : choices[0])?.focus();
    setPendingTab(null);
  }, [pendingTab, callOpen]);
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
    if (sampleId === undefined) return;
    history.record(customIndexEntry(ops, {
      label: `added a custom ${c.phase} index on ${sample?.name ?? "this sample"}`, sampleId, exposureId: exposure.id,
      phase: c.phase, basis: c.basis, ratios: c.ratios, wasFormFactor: isFormFactor,
    }));
  }, [draft, exposure, canCall, commitCustom, sampleId, sample?.name, isFormFactor, history, ops]);
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
      customPanel.current?.focus();
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
    const hovered = candidates.filter((c) => hoverIds?.includes(c.id));
    const shown = hovered.length > 0 ? hovered : candidates.filter((c) => assigned.has(c.id));
    return shown.map((c) => {
      const fit = storedFit(c, xExtent);
      return { key: String(c.id), phase: c.phase, color: phaseColor(c.phase), fit, lattice: c.lattice_d, param: latticeParam(c.phase), r2: c.r_squared ?? fit.r2 };
    });
  }, [draft, draftFit, xExtent, candidates, hoverIds, assigned]);

  const applyHistory = useCallback(
    (dir: "undo" | "redo") => {
      if (!username) return;
      if (draft) {
        showToast("Close the custom index first (Esc).");
        return;
      }
      const saving = ops.saving();
      const step: UndoStep | undefined = dir === "undo" ? history.undo() : history.redo();
      if (!step) {
        showToast(dir === "undo" ? "Nothing to undo." : "Nothing to redo.");
        return;
      }
      const { entry } = step;
      if (entry.sampleId !== sampleId) {
        const qs = entry.exposureId === undefined ? "" : `?exposure=${entry.exposureId}`;
        navigate(`/experiments/${experimentId}/workspace/${entry.sampleId}${qs}`, { replace: true });
      } else if (entry.exposureId !== undefined) {
        setExposureId(entry.exposureId);
      }
      const verb = dir === "undo" ? "Undone" : "Redone";
      if (saving) showToast(`Saving your last edit first, then ${dir === "undo" ? "undoing" : "redoing"}: ${entry.label}`);
      step.done.then(
        (sent) => showToast(
          sent > 0 ? `${verb}: ${entry.label}` : `Already that way, nothing to ${dir}: ${entry.label}`,
        ),
        () => showToast(`Couldn't ${dir}: ${entry.label}. The server can't be reached.`, "error"),
      );
    },
    [username, draft, ops, history, sampleId, navigate, experimentId],
  );

  useEffect(() => {
    const onKey = (ev: KeyboardEvent): void => {
      if ((ev.metaKey || ev.ctrlKey) && !ev.altKey && !ev.shiftKey && ev.key === "0" && !suppressGlobalKeys(ev)) {
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
      const lower = ev.key.toLowerCase();
      const redoY = ev.ctrlKey && !ev.metaKey && !ev.shiftKey && lower === "y";
      if ((ev.metaKey || ev.ctrlKey) && !ev.altKey && (lower === "z" || redoY) && !suppressGlobalKeys(ev)) {
        ev.preventDefault();
        if (ev.repeat) return;
        applyHistory(redoY || ev.shiftKey ? "redo" : "undo");
        return;
      }
      if (helpKeyFromEmptyField(ev)) {
        ev.preventDefault();
        setHelpOpen(true);
        return;
      }
      if (ev.metaKey || ev.ctrlKey || ev.altKey) return;
      if (suppressGlobalKeys(ev)) {
        if (ev.key === "Escape" && !ev.defaultPrevented && !inDialog(ev) && isEditable(ev.target)) (ev.target as HTMLElement).blur();
        return;
      }
      if (ev.key === "/") {
        ev.preventDefault();
        const filter = document.getElementById("workspace-filter");
        filter?.focus();
        if (filter instanceof HTMLInputElement) filter.select();
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
        } else if (ev.key === "Tab") {
          ev.preventDefault();
          setDraftClamped({ ...draft, phase: stepSymmetry(draft.phase, ev.shiftKey ? -1 : 1) });
        } else if (ev.key === "?") {
          ev.preventDefault();
          setHelpOpen(true);
        } else if (ev.key === " " && !ev.repeat && tile && imageUrl(tile, false)) {
          if (ev.target instanceof HTMLElement && ev.target.closest("button,a")) return;
          ev.preventDefault();
          setBigOpen(true);
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
        case "Tab": {
          if (tabFree.current) {
            handled = false;
            break;
          }
          const next = cycleChoice(callChoices(), ev.target, ev.shiftKey);
          if (next) next.focus();
          else if (!callOpen) {
            setOpenPanel("call");
            setPendingTab(ev.shiftKey ? "last" : "first");
          } else handled = false;
          break;
        }
        case "Escape":
          if (inCallChoices(ev.target)) {
            tabFree.current = true;
            (ev.target as HTMLElement).blur();
          } else handled = false;
          break;
        case "ArrowDown":
          step(1);
          break;
        case "ArrowUp":
          step(-1);
          break;
        case "ArrowRight":
        case "ArrowLeft": {
          const eid = stepExposure(exposureList, tile?.id, ev.key === "ArrowRight" ? 1 : -1);
          if (eid !== undefined) setExposureId(eid);
          break;
        }
        case "[":
          setSampleListHidden(!sampleListHidden);
          break;
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
        case "v":
        case "V":
          setView(stepView(view, ev.shiftKey ? -1 : 1));
          break;
        case "z":
        case "Z": {
          const nid = nextUnreviewed(visible, sampleId, 1);
          if (nid !== undefined) goSample(nid);
          break;
        }
        case " ":
          if (ev.repeat) break;
          if (tile && imageUrl(tile, false)) setBigOpen(true);
          else handled = false;
          break;
        case "x":
        case "X":
          if (ev.repeat) break;
          if (tile && username) toggleDrop(tile);
          break;
        case "Enter": {
          const onControl = ev.target instanceof HTMLElement
            && ev.target.closest("button,a,[role=checkbox]") !== null
            && ev.target.closest('[data-roving="image"]') === null;
          if (onControl || !tileSelectable) handled = false;
          else if (!ev.repeat) toggleTileSelection();
          break;
        }
        case "0":
          if (username) {
            setOpenPanel("call");
            toggleFormFactor();
          }
          break;
        default: {
          const n = Number(ev.key);
          const row = Number.isInteger(n) && n >= 1 ? keyed[n - 1] : undefined;
          if (row && username) {
            setOpenPanel("call");
            if (row.kind === "group") toggleGroup(row.members);
            else toggleCandidate(row.candidate);
          } else handled = false;
        }
      }
      if (handled) ev.preventDefault();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [visible, sampleId, goSample, exposureList, tile, toggleDrop, tileSelectable, toggleTileSelection, keyed, toggleCandidate, toggleGroup, username, draft, addCustom, setDraftClamped, openCustom, toggleFormFactor, view, setView, sample, openLabeling, applyHistory, callOpen, sampleListHidden, setSampleListHidden]);

  const droppedCount = ownFrames.filter((e) => e.status === "rejected").length;
  const keptHint = droppedCount > 0 ? `${ownFrames.length - droppedCount} of ${ownFrames.length} kept` : undefined;
  const reviewed = samples.filter((s) => sampleStatus(s).kind !== "unreviewed").length;
  const countOf = (at: typeof tileFrame): string | undefined => (at ? `Image ${at.index + 1} of ${at.group.exposures.length}` : undefined);
  const tileOwnerSample = tileFrame ? inExperiment.find((s) => s.id === tileFrame.group.sampleId) : undefined;
  const tileOwner = tileOwnerSample ? frameSetupName(tileOwnerSample) : undefined;
  const tileCount = countOf(tileFrame);
  const frameLine = tileCount && tileFrame && tileFrame.group.sampleId !== sampleId ? `${tileCount}, ${tileOwner ?? ""}` : tileCount;
  const emptyMessage = !Number.isInteger(experimentId)
    ? "This workspace link is not valid."
    : corpus.isLoading || !setupsKnown
      ? "Loading samples…"
      : samples.length === 0
        ? ingesting ? "No samples yet. They appear as the ingest reads the data folder." : "No samples in this experiment."
        : inExperiment.some((s) => s.id === requested)
          ? "This sample's setup is hidden. Show it again from Settings in the palette."
          : "This sample is not in this experiment.";

  const addReference = (rid: number): void => {
    if (sampleId === undefined) return;
    if (comparesFull) {
      showToast(COMPARES_FULL);
      return;
    }
    patchCompares.mutate({ sampleId, change: { add: [rid] } });
    history.record(compareEntry(ops, {
      label: `compared ${sample?.name ?? "this sample"} with ${nameOf(rid)}`, sampleId, exposureId: undefined,
      edit: { kind: "add", ref: rid },
    }));
    recordRecent(recentKind, String(rid));
    setRecentCompares(recentRank(recentKind));
    setShown(rid, true);
    setRefsOn(true);
  };
  const toggleSubtract = (rid: number): void => {
    if (sampleId === undefined) return;
    setRefsOn(true);
    setMinus(subRef !== rid);
    const to = subRef === rid ? undefined : rid;
    patchCompares.mutate({ sampleId, change: { subtraction: to ?? null } });
    history.record(compareEntry(ops, {
      label: to === undefined ? `stopped subtracting ${nameOf(rid)} from ${sample?.name ?? "this sample"}` : `subtracted ${nameOf(rid)} from ${sample?.name ?? "this sample"}`,
      sampleId, exposureId: undefined, edit: { kind: "subtract", from: subChoice, to },
    }));
  };
  const toggleCompare = (rid: number): void => {
    if (rid === sampleId) return;
    if (refs.includes(rid)) removeCompare(rid);
    else addReference(rid);
  };
  const palette = workspacePalette({
    setupOf: (id) => sampleSetup.get(id),
    modKey: MOD_KEY,
    canEdit: Boolean(username && sample),
    nextUnreviewed: () => {
      const nid = nextUnreviewed(visible, sampleId, 1);
      if (nid !== undefined) goSample(nid);
    },
    compare: {
      samples: pickable,
      notes: bufferNotes,
      pick: addReference,
      full: comparesFull ? COMPARES_FULL : undefined,
      compared: refs.filter((r) => r !== sampleId).flatMap((r) => samples.filter((x) => x.id === r)),
      remove: removeCompare,
      subtractable: activeRefs.flatMap((r) => samples.filter((x) => x.id === r)),
      subtracted: subRef,
      subtract: toggleSubtract,
    },
    bridge: {
      samples: bridgePickable,
      pick: (rid, keep) => {
        if (sampleId === undefined) return;
        bridgeTo(rid, keep);
      },
      choices: bridgeChoices,
      names: (id) => nameWithSetup({ id, name: sampleNames.get(id) ?? `Sample ${id}` }),
      bridged: bridgedNames,
      leave: leaveBridges,
      groupChoices: (ids) => (sample ? groupKeepChoices(sample, ids.flatMap((id) => samples.filter((s) => s.id === id)), nameWithSetup) : undefined),
      pickMany: (ids, keep) => void bridgeMany(ids, keep),
    },
    candidates: {
      list: choices,
      ticked: assigned,
      tick: (ix) => {
        setOpenPanel("call");
        toggleCandidate(ix);
      },
    },
    labels: {
      known: paletteLabels,
      colorFor: (name) => (labelColorIsDefault(labelColors.data ?? [], name) ? firstUnusedColor(labelColors.data ?? [], vocabulary) : labelColor(labelColors.data ?? [], name)),
      toggle: (l) => editLabel(l.label, l.onSample !== "confirmed"),
      create: createLabel,
      recolour: recolorLabel,
    },
    customIndex: openCustom,
    ...(tileSelectable ? { image: { selected: tileSelected, toggle: toggleTileSelection } } : {}),
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
    ...(difference && subRef !== undefined ? { difference: { on: minusOn, name: subName, toggle: () => setMinus(!minusOn) } } : {}),
    settings:
      username && experiment.data
        ? {
            name: experiment.data.name ?? "",
            description: experiment.data.description ?? "",
            save: editExperiment,
            setups: paletteSetups,
            saveSetup: editSetup,
          }
        : undefined,
    rescan: ingesting ? undefined : () => triggerScan.mutate(ingestFailed || failedAnalyses > 0),
    problems: problemRows.length > 0 || unrecognizedSetups > 0 || ingestFailed ? () => setProblemsOpen(true) : undefined,
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
    history: {
      undo: username ? history.nextUndo?.label : undefined,
      redo: username ? history.nextRedo?.label : undefined,
      run: applyHistory,
    },
    samples,
    goSample,
    elsewhere: { ...elsewhere, go: (path) => navigate(path) },
  });

  const customIndexButton = <CustomIndexButton open={draft !== null} disabled={!canCall} onClick={draft ? () => setDraft(null) : openCustom} />;
  const indexingSection = (
    <Section
      title="Indexing"
      className="shrink-0"
    >
      <IndexingPlot series={indexSeries} />
    </Section>
  );

  return (
    <WorkspaceLayout
      theme={theme.theme}
      navHidden={sampleListHidden}
      onShowNav={() => setSampleListHidden(false)}
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
              onAddSetups={() => navigate(`/experiments/new?path=${encodeURIComponent(experiment.data?.path ?? "")}`)}
              onClose={() => setProblemsOpen(false)}
            />
          )}
          {bigOpen && tile && (
            <DetectorOverlay
              exposures={exposureList}
              currentId={tile.id}
              caption={frameLine}
              onStep={(dir) => {
                const eid = stepExposure(exposureList, tile.id, dir);
                if (eid !== undefined) setExposureId(eid);
              }}
              onToggleDrop={username ? () => toggleDrop(tile) : undefined}
              selection={tileSelectable ? { selected: tileSelected, toggle: toggleTileSelection } : undefined}
              onClose={() => setBigOpen(false)}
              ringsOf={ringsOf}
            />
          )}
          {paletteOpen && <Palette groups={palette.groups} pages={palette.pages} start={paletteStart} onClose={() => { setPaletteOpen(false); setPaletteStart(undefined); }} onType={paletteCorpus.onType} />}
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
      nav={<SampleList samples={visible} setupHeaders={setupHeaders} currentId={sampleId} filter={filter} onFilter={setFilter} onPick={goSample} onToggleCompare={toggleCompare} compareColors={compareColors} hiddenCompares={shownCompareSet} bridgeCounts={counts} listHidden={sampleListHidden} onToggleList={() => setSampleListHidden(!sampleListHidden)} notice={<ExperimentNotice progress={ingestProgress} ingestFailed={ingestFailed} waiting={waitingForTrace} failed={failedAnalyses} unrecognized={unrecognizedSetups} onDetails={() => setProblemsOpen(true)} />} />}
      main={
        sample ? (
          <>
            <SampleHeading
              name={sample.name}
              badge={bridgeMates.length > 0 ? <BridgeMenu key={sampleId} mates={bridgeMates} editable={Boolean(username)} onGo={goSample} onRemove={unbridgeOne} /> : undefined}
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
                  onSetColor={recolorLabel}
                  onAdd={(label) => editLabel(label, true)}
                  onRemove={(label) => editLabel(label, false)}
                />
              }
            />
            <div className="relative px-3.5 py-1">
              <CompareBar
                chips={activeRefs.map((rid) => ({
                  id: rid,
                  name: sampleNames.get(rid) ?? `Sample ${rid}`,
                  setup: sampleSetup.get(rid) !== (sampleId === undefined ? undefined : sampleSetup.get(sampleId)) ? sampleSetup.get(rid) : undefined,
                  color: COMPARE_COLORS[activeRefs.indexOf(rid) % COMPARE_COLORS.length]!,
                  scale: scales.get(rid),
                  hidden: shownCompareSet.has(rid),
                }))}
                subtractId={subRef}
                hidden={!refsOn}
                modKey={MOD_KEY}
                onToggleShown={(rid) => {
                  if (minusOn && rid !== subRef) {
                    setMinus(false);
                    setShown(rid, true);
                    return;
                  }
                  setRefsOn(true);
                  setShown(rid, !(refsOn && !hiddenRefs.has(rid)));
                }}
                onSubtract={toggleSubtract}
                onRemove={removeCompare}
                onAdd={() => setPicking(true)}
                full={comparesFull}
                trailing={difference && subRef !== undefined ? <SubtractToggle name={subName} on={minusOn} onChange={setMinus} /> : undefined}
              />
              {picking && (
                <ComparePicker
                  samples={pickable}
                  notes={bufferNotes}
                  setups={sampleSetup}
                  full={comparesFull ? COMPARES_FULL : undefined}
                  onPick={addReference}
                  onClose={() => setPicking(false)}
                />
              )}
            </div>
            <div className="mx-3.5">
              <PlotCard {...(darkNotice ? { notice: darkNotice } : {})} corner={<PlotViewToggle options={PLOT_VIEWS} view={view} onChange={setView} />}>
                <WorkspacePlot
                  trace={trace.data}
                  {...(noTraceMessage
                    ? { emptyMessage: noTraceMessage }
                    : trace.isError
                      ? { onRetry: () => void trace.refetch(), retrying: trace.isFetching }
                      : {})}
                  peaks={peaks.data ?? []}
                  combs={combs}
                  highlight={highlight}
                  compares={compares}
                  members={members}
                  {...(minusOn && difference ? { difference } : {})}
                  xDomain={shownDomain}
                  view={view}
                  onXDomain={setXDomain}
                  {...(draft
                    ? { onAddPeak: snapDraft, onClickPeak: snapDraftToPeak }
                    : canEditPeaks
                      ? { onAddPeak: addPeakAt, onClickPeak }
                      : {})}
                />
              </PlotCard>
            </div>
            <MainSection title="Images" hint={keptHint}>
              <ExposureStrip groups={groups.map((g) => ({ key: g.sampleId, title: g.title, exposures: g.exposures, selection: username ? groupSelections.get(g.sampleId) : undefined }))} currentId={tile?.id} onPick={setExposureId} onToggleDrop={toggleDrop} warningOf={beamOffWarning} />
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
          <div className={callOpen ? "flex min-h-24 flex-1 flex-col overflow-y-auto" : "shrink-0"} data-testid="workspace-tools-middle">
            {draft && draftBounds ? (
              <Section title="Assignment" aside={customIndexButton} pinAside collapsible={shortWindow ? { open: true } : undefined}>
                <CustomIndexPanel
                  panelRef={customPanel}
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
              <Section title="Assignment" aside={customIndexButton} pinAside collapsible={shortWindow ? { open: callOpen, onOpen: () => setOpenPanel("call") } : undefined}>
                <CandidateList
                  rows={rows}
                  assigned={assigned}
                  unsupported={unsupportedIds}
                  disabled={!username}
                  formFactor={isFormFactor}
                  onToggle={toggleCandidate}
                  onToggleGroup={toggleGroup}
                  onFormFactor={toggleFormFactor}
                  onHover={hoverCandidates}
                />
              </Section>
            )}
          </div>
          {!shortWindow && indexingSection}
          <Section
            title="Image"
            aside="← →"
            folded={groups.length > 1 ? tileOwner : undefined}
            className={shortWindow && exposureOpen ? "min-h-0 flex-1 overflow-y-auto" : exposureOpen ? "min-h-0 shrink overflow-y-auto" : "shrink-0"}
            collapsible={shortWindow ? { open: exposureOpen, onOpen: draft === null ? () => setOpenPanel("exposure") : undefined } : undefined}
          >
            <MetaList rows={exposureRows(tile, tileOwner)} />
          </Section>
          {shortWindow && indexingSection}
        </>
      }
    />
  );
}
