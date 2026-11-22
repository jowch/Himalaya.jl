import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, LabelColor, SampleCallMember } from "../api";
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
} from "../queries";
import { comparesOf, EMPTY_COMPARES } from "../lib/compares";
import { byRecent, recentRank, recordRecent } from "../lib/recents";
import { helpKeyFromEmptyField, inDialog, isEditable, suppressGlobalKeys } from "../lib/keys";
import { useUndoHistory, type UndoEntry, type UndoStep } from "./undo";
import {
  bridgeEntry, combinedEntry, settleCallEntry, customIndexEntry, dropEntry, excludeEntry, identify, labelEntry, noPhaseEntry, peakEntry,
  experimentEntry, setupEntry, labelColorEntry, tickEntry, unbridgeEntry, untickEntry, compareEntry,
} from "./undoEntries";
import { showToast } from "../lib/toast";
import { useScanState } from "../hooks/useScanState";
import { useMediaQuery } from "./useMediaQuery";
import { isIngestStage, stageLabel } from "../lib/ingestStages";
import { bridgeCounts, bridgedWith } from "../lib/bridges";
import { callStatesSplit, keepChoices, type KeepCall } from "../lib/bridgeCall";
import { bufferRole, buffersFirst, displayLabel, firstUnusedColor, labelColor, labelColorIsDefault, labelsFor, sameLabel, sampleMatchesFilter, labelVocabulary, labelKey } from "../lib/labels";
import { representativeExposure } from "../lib/exposure";
import { useCorpusPalette } from "../hooks/useCorpusPalette";
import { allFrames, frameGroups, frameOf } from "./frameGroups";
import { phaseColor } from "../phases";
import { useAppState } from "../state";
import { firstUnreviewed, i0BelowZero, traceIsDark, imageUrl, nearestKeptExposure, nextUnreviewed, peakClick, rankCandidates, sampleStatus, stepExposure } from "./model";
import { samePeak, type PeakRef } from "../lib/peakRef";
import { clampDomain, compareScale, COMPARES_FULL, MAX_COMPARES, scaledTrace, subtractTrace } from "./compare";
import type { CompareScale } from "./compare";
import { useCompareTraces } from "./useCompareTraces";
import { usePrefetchNeighbours } from "./prefetch";
import { CompareBar } from "../ui/CompareBar";
import { ComparePicker } from "../ui/ComparePicker";
import { SubtractToggle } from "../ui/SubtractToggle";
import { PlotCard } from "../ui/PlotCard";
import { HelpOverlay, type HelpSection } from "../ui/HelpOverlay";
import { MOD_KEY } from "../lib/modKey";
import { DetectorOverlay } from "../ui/DetectorOverlay";
import { Palette } from "../ui/palette/Palette";
import { keepCallForm, workspacePalette, type PaletteSetup } from "./palette";
import type { PaletteForm } from "../ui/palette/model";
import { useAppTheme } from "../lib/theme";
import { PLOT_VIEWS, usePlotView } from "./plotView";
import { positiveExtent } from "../plot/projection";
import { CustomIndexPanel, NUDGE_STEP, NUDGE_STEP_LARGE } from "../ui/CustomIndexPanel";
import { IndexingPlot, type IndexSeries } from "../ui/IndexingPlot";
import { clampLattice, draftCommit, draftPrediction, indexFit, initialDraft, latticeParam, storedFit, type CustomDraft } from "./indexing";
import { latticeBounds, latticeForFirstOrderOnPeak } from "../lib/customIndex";

import { WorkspaceHeader } from "../ui/WorkspaceHeader";
import { Note } from "../ui/Note";
import { WorkspaceLayout } from "../ui/WorkspaceLayout";
import { SampleList } from "../ui/SampleList";
import { ExperimentNotice } from "../ui/ExperimentNotice";
import { ProblemList } from "../ui/ProblemList";
import { RemoveExperimentDialog } from "../ui/RemoveExperimentDialog";
import { COMPARE_COLORS, WorkspacePlot, type CombLine, type CompareLine } from "../ui/WorkspacePlot";
import { ExposureStrip } from "../ui/ExposureStrip";
import { CandidateList } from "../ui/CandidateList";
import { callRows, keyedRows } from "./bonnetGroups";
import { MetaList } from "../ui/MetaList";
import { MainSection, SampleHeading, Section } from "../ui/Section";
import { claimedPeaks } from "./claimedPeaks";
import { CustomIndexButton } from "../ui/CustomIndexButton";
import { LabelRow } from "../ui/LabelRow";
import type { ExperimentPatch } from "../lib/queue/mutators/experiment";
import type { SetupPatch } from "../lib/queue/mutators/setup";

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
      title: "Frames",
      rows: [
        { keys: ["←", "→"], what: "Previous or next frame" },
        { keys: ["X"], what: "Drop the frame and show the nearest kept one, or restore a dropped one" },
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

const EXPOSURE_LABELS = ["Sample", "File", "Kind", "Acquired", "Exposure", "Flight tube", "Energy", "Detector", "I₀", "I signal"];

function exposureRows(e: Exposure | undefined, owner: string | undefined): [string, string][] {
  if (!e) return EXPOSURE_LABELS.map((label) => [label, "—"]);
  const rows: [string, string][] = [["Sample", owner ?? "—"], ["File", e.filename ?? "—"]];
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
  const bridged = useMemo(() => (sampleId === undefined ? [] : bridgedWith(bridgeList, sampleId)), [bridgeList, sampleId]);
  const exposures = useExposures(sampleId);
  const mateExposures = useSamplesExposures(bridged);
  const [exposureId, setExposureId] = useState<number | undefined>(linkedExposure);
  useEffect(() => setExposureId(linkedExposure), [sampleId, linkedExposure]);
  const ownFrames: Exposure[] = useMemo(() => exposures.data ?? [], [exposures.data]);
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
      return s ? [{ sample: s, frames: mateExposures[i] ?? [] }] : [];
    });
    return frameGroups({ sample, frames: ownFrames }, mates, (s) => setupOrder.get(s.setup_id ?? -1) ?? 0, frameSetupName);
  }, [sample, bridged, inExperiment, ownFrames, mateExposures, setupOrder, frameSetupName]);
  const exposureList: Exposure[] = useMemo(() => allFrames(groups), [groups]);
  const exposure = exposureList.find((e) => e.id === exposureId) ?? representativeExposure(ownFrames);
  const shownFrame = frameOf(groups, exposure?.id);

  const trace = useTrace(exposure?.id);
  const noTrace = trace.error instanceof ApiError && trace.error.status === 404;
  const i0Hint = i0BelowZero(exposure);
  const darkNotice = traceIsDark(trace.data)
    ? `No counts on the detector: the beam was likely off or the shutter closed.${i0Hint ? " I₀ also reads below zero." : ""}`
    : undefined;
  const tracePending = problemRows.some((p) => p.problem === "no_trace" && p.exposure_id === exposure?.id);
  const noTraceMessage = !noTrace
    ? undefined
    : tracePending
      ? "Waiting for this frame's trace."
      : `No integrated trace for this frame.${i0Hint ? " I₀ reads below zero, so the beam was likely off." : ""}`;
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
  const rows = useMemo(() => callRows(candidates), [candidates]);
  const keyed = useMemo(() => keyedRows(rows), [rows]);
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
    for (const c of candidates) if (phases.has(c.phase)) here.add(c.id);
    return here;
  }, [assignment.data, otherMembers, candidates]);
  const [hoverIds, setHoverIds] = useState<readonly number[] | undefined>(undefined);

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
  const ops = useEditOps(experimentId);
  const history = useUndoHistory(`${experimentId}:${username ?? ""}`, ops.whenSaved);
  const frameName = useCallback(
    (eid: number) => {
      const at = frameOf(groups, eid);
      const owner = at && at.group.sampleId !== sample?.id ? inExperiment.find((s) => s.id === at.group.sampleId)?.name : sample?.name;
      return `frame ${(at?.index ?? -1) + 1} of ${owner ?? "this sample"}`;
    },
    [groups, sample?.name, sample?.id, inExperiment],
  );

  const [xDomain, setXDomain] = useState<[number, number] | null>(null);
  const [view, setView] = usePlotView();
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
  useEffect(() => setDraft(null), [exposure?.id]);
  const commitCustom = useCommitCustomIndex(exposure?.id ?? 0);
  const bridgedNames = useMemo(() => {
    const names = new Map(inExperiment.map((s) => [s.id, s.name]));
    return bridged.map((r) => names.get(r) ?? `Sample ${r}`);
  }, [bridged, inExperiment]);
  const bridgedTitles = useMemo(
    () => bridged.map((r) => {
      const s = inExperiment.find((x) => x.id === r);
      return s ? frameSetupName(s) : `Sample ${r}`;
    }),
    [bridged, inExperiment, frameSetupName],
  );
  const counts = useMemo(() => bridgeCounts(bridgeList), [bridgeList]);
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

  const onClickPeak = useCallback(
    (ref: PeakRef) => {
      const action = peakClick(peaks.data ?? [], ref);
      if (!action || !exposure || sampleId === undefined) return;
      const eid = exposure.id;
      const q = (peaks.data ?? []).find((pk) => samePeak(pk, ref))?.q;
      if (action.kind === "remove") {
        if (ref.id < 0) return;
        removePeak.mutate(ref.id);
        if (q !== undefined) {
          history.record(peakEntry(ops, {
            label: `removed a peak at q = ${q.toFixed(4)} on ${frameName(eid)}`, sampleId, exposureId: eid, q, added: false,
          }));
        }
        return;
      }
      setExcluded.mutate({ peakId: ref.id, excluded: action.excluded });
      if (q !== undefined) {
        history.record(excludeEntry(ops, {
          label: `${action.excluded ? "turned off" : "turned on"} a peak on ${frameName(eid)}`, sampleId, exposureId: eid, q,
          excluded: action.excluded,
        }));
      }
    },
    [peaks.data, removePeak, setExcluded, exposure, sampleId, history, ops, frameName],
  );

  const addPeakAt = useCallback(
    (q: number) => {
      if (!exposure || sampleId === undefined) return;
      const eid = exposure.id;
      addPeak.mutate(q);
      history.record(peakEntry(ops, {
        label: `added a peak at q = ${q.toFixed(4)} on ${frameName(eid)}`, sampleId, exposureId: eid, q, added: true,
      }));
    },
    [exposure, sampleId, addPeak, history, ops, frameName],
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
      else if (e.id === exposure?.id) setExposureId(nearestKeptExposure(exposureList, e.id) ?? e.id);
      if (sampleId === undefined) return;
      const owner = e.sample_id;
      if (owner === sampleId) setStatus.mutate({ exposureId: e.id, status: next });
      else ops.setStatus(owner, e.id, next);
      history.record(dropEntry(ops, {
        label: `${next ? "dropped" : "restored"} ${frameName(e.id)}`, sampleId, exposureId: e.id, owner, prev, next,
      }));
    },
    [setStatus, sampleId, history, ops, frameName, exposureList, exposure],
  );

  const draftQs = useMemo(() => (draft ? draftPrediction(draft) : []), [draft]);
  const combs: CombLine[] = useMemo(() => {
    if (draft) return [{ id: -1, qs: draftQs, color: phaseColor(draft.phase) }];
    return candidates.filter((c) => hoverIds?.includes(c.id)).map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase) }));
  }, [candidates, hoverIds, draft, draftQs]);

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
    compareTraces.forEach((t, rid) => out.set(rid, compareScale(mainTrace, t, peakQs)));
    return out;
  }, [mainTrace, compareTraces, peakQs]);
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
  const bridgePickable = useMemo(() => samples.filter((s) => s.id !== sampleId && !bridged.includes(s.id)), [samples, sampleId, bridged]);
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
  const setupNames = useMemo(() => (setups.length > 1 ? new Map(setups.map((st) => [st.id, st.name])) : undefined), [setups]);
  const setupOfSample = useCallback(
    (s: CorpusSample | undefined): string | undefined => (s?.setup_id == null ? undefined : setupNames?.get(s.setup_id)),
    [setupNames],
  );
  const nameWithSetup = useCallback(
    (s: Pick<CorpusSample, "id" | "name">): string => {
      const setup = setupOfSample(samples.find((x) => x.id === s.id));
      return setup ? `${s.name} in ${setup}` : s.name;
    },
    [samples, setupOfSample],
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
          if (ev.repeat) break;
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
  }, [visible, sampleId, goSample, exposureList, exposure, toggleDrop, keyed, toggleCandidate, toggleGroup, username, draft, addCustom, setDraftClamped, openCustom, toggleFormFactor, view, setView, sample, openLabeling, applyHistory]);

  const droppedCount = ownFrames.filter((e) => e.status === "rejected").length;
  const keptHint = droppedCount > 0 ? `${ownFrames.length - droppedCount} of ${ownFrames.length} kept` : undefined;
  const reviewed = samples.filter((s) => sampleStatus(s).kind !== "unreviewed").length;
  const ownerSample = shownFrame ? inExperiment.find((s) => s.id === shownFrame.group.sampleId) : undefined;
  const frameOwner = ownerSample ? frameSetupName(ownerSample) : undefined;
  const frameCount = shownFrame ? `Frame ${shownFrame.index + 1} of ${shownFrame.group.exposures.length}` : undefined;
  const frameLine = frameCount && shownFrame && shownFrame.group.sampleId !== sampleId ? `${frameCount}, ${frameOwner ?? ""}` : frameCount;
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
      setupOf: (id) => setupOfSample(samples.find((x) => x.id === id)),
      bridged: bridgedNames,
      leave: leaveBridges,
    },
    candidates: {
      list: candidates,
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

  const indexingSection = (
    <Section
      title="Indexing"
      className="shrink-0"
      aside={<CustomIndexButton open={draft !== null} disabled={!canCall} onClick={draft ? () => setDraft(null) : openCustom} />}
    >
      <IndexingPlot series={indexSeries} />
    </Section>
  );

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
              onAddSetups={() => navigate(`/experiments/new?path=${encodeURIComponent(experiment.data?.path ?? "")}`)}
              onClose={() => setProblemsOpen(false)}
            />
          )}
          {bigOpen && exposure && (
            <DetectorOverlay
              exposures={exposureList}
              currentId={exposure.id}
              caption={frameLine}
              onStep={(dir) => {
                const eid = stepExposure(exposureList, exposure.id, dir);
                if (eid !== undefined) setExposureId(eid);
              }}
              onToggleDrop={username ? () => toggleDrop(exposure) : undefined}
              onClose={() => setBigOpen(false)}
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
      nav={<SampleList samples={visible} setupHeaders={setupHeaders} currentId={sampleId} filter={filter} onFilter={setFilter} onPick={goSample} onToggleCompare={toggleCompare} compareColors={compareColors} hiddenCompares={shownCompareSet} bridgeCounts={counts} notice={<ExperimentNotice progress={ingestProgress} ingestFailed={ingestFailed} waiting={waitingForTrace} failed={failedAnalyses} unrecognized={unrecognizedSetups} onDetails={() => setProblemsOpen(true)} />} />}
      main={
        sample ? (
          <>
            <SampleHeading
              name={sample.name}
              meta={exposure ? frameCount : undefined}
              aside={
                <LabelRow
                  note={bridgedTitles.length > 0 ? `Same as ${bridgedTitles.join(", ")}` : undefined}
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
                trailing={difference && subRef !== undefined ? <SubtractToggle name={subName} on={minusOn} onChange={setMinus} /> : undefined}
              />
              {picking && (
                <ComparePicker
                  samples={pickable}
                  notes={bufferNotes}
                  full={comparesFull ? COMPARES_FULL : undefined}
                  onPick={addReference}
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
                  highlight={highlight}
                  compares={compares}
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
            <MainSection title="Frames" hint={keptHint}>
              <ExposureStrip groups={groups.map((g) => ({ key: g.sampleId, title: g.title, exposures: g.exposures }))} currentId={exposure?.id} onPick={setExposureId} onToggleDrop={toggleDrop} />
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
                  rows={rows}
                  assigned={assigned}
                  disabled={exposureDropped || !username}
                  formFactor={isFormFactor}
                  onToggle={toggleCandidate}
                  onToggleGroup={toggleGroup}
                  onFormFactor={toggleFormFactor}
                  onHover={setHoverIds}
                />
              </Section>
            )}
          </div>
          {!shortWindow && indexingSection}
          <Section
            title="Frame"
            aside="← →"
            className={shortWindow && exposureOpen ? "min-h-0 flex-1 overflow-y-auto" : exposureOpen ? "min-h-0 shrink overflow-y-auto" : "shrink-0"}
            collapsible={shortWindow ? { open: exposureOpen, onOpen: draft === null ? () => setOpenPanel("exposure") : undefined } : undefined}
          >
            <MetaList rows={exposureRows(exposure, frameOwner)} />
          </Section>
          {shortWindow && indexingSection}
        </>
      }
    />
  );
}
