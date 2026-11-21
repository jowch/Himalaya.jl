import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import type { Exposure, IndexEntry } from "../api";
import {
  useAddAssignmentPhase,
  useAddPeak,
  useCommitCustomIndex,
  useAssignment,
  useCorpusSamples,
  useExperiment,
  useExposures,
  useIndices,
  usePeaks,
  useRemoveAssignmentPhase,
  useSetAssignmentState,
  useRemovePeak,
  useSelectExposure,
  useSetExposureStatus,
  useSetPeakExcluded,
  useTrace,
} from "../queries";
import { representativeExposure } from "../lib/exposure";
import { phaseColor } from "../phases";
import { useAppState } from "../state";
import { nextUnreviewed, peakClick, rankCandidates, sampleStatus, stepExposure } from "./model";
import { addCompare, clampDomain, compareScale, scaledTrace, subtractTrace } from "./compare";
import type { CompareScale } from "./compare";
import { useCompareTraces } from "./useCompareTraces";
import { CompareBar } from "./ui/CompareBar";
import { ComparePicker } from "./ui/ComparePicker";
import { SubtractionStrip } from "./ui/SubtractionStrip";
import { HelpOverlay, type HelpSection } from "./ui/HelpOverlay";
import { CommandPalette, type Command } from "./ui/CommandPalette";
import { useWorkspaceTheme } from "./theme";
import { PLOT_VIEWS, nextView, usePlotView } from "./plotView";
import { positiveExtent } from "../print/plot/projection";
import { CustomIndexPanel, NUDGE_STEP, NUDGE_STEP_LARGE } from "./ui/CustomIndexPanel";
import { IndexingPlot, type IndexSeries } from "./ui/IndexingPlot";
import { clampLattice, draftCommit, draftPrediction, indexFit, initialDraft, storedFit, type CustomDraft } from "./indexing";
import { latticeBounds, latticeForFirstOrderOnPeak } from "../lib/customIndex";

import { WorkspaceHeader } from "./ui/WorkspaceHeader";
import { Note } from "./ui/Note";
import { WorkspaceLayout } from "./ui/WorkspaceLayout";
import { SampleList } from "./ui/SampleList";
import { COMPARE_DASHES, WorkspacePlot, type CombLine, type CompareLine } from "./ui/WorkspacePlot";
import { ExposureStrip } from "./ui/ExposureStrip";
import { CandidateList } from "./ui/CandidateList";
import { MetaList } from "./ui/MetaList";
import { Section } from "./ui/Section";

function helpSections(mod: string): HelpSection[] {
  return [
    {
      title: "Move",
      rows: [
        { keys: ["↑", "↓"], what: "Previous or next sample in the list" },
        { keys: ["J", "Shift J"], what: "Next or previous unreviewed sample" },
        { keys: ["/"], what: "Filter samples by name" },
      ],
    },
    {
      title: "Exposures",
      rows: [
        { keys: ["←", "→"], what: "Previous or next exposure" },
        { keys: ["X"], what: "Drop or restore the exposure" },
      ],
    },
    {
      title: "Call",
      rows: [
        { keys: ["1–9"], what: "Tick or untick a phase candidate" },
        { keys: ["0"], what: "Mark no phase (form factor), or clear it" },
        { keys: ["I"], what: "Custom index: pick a symmetry and drag the lattice; Enter adds, Esc cancels" },
      ],
    },
    {
      title: "Plot",
      rows: [
        { gesture: "Click", what: "Turn a peak off or on, remove one you added, or add one" },
        { gesture: "Scroll", what: "Zoom q" },
        { keys: ["Z"], what: "Reset the zoom" },
        { keys: ["L", "Shift L"], what: "Next or previous view: log, linear, log-log, q²I" },
        { keys: ["C"], what: "Compare with another sample" },
        { keys: ["R"], what: "Hide or show comparisons" },
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
  const experiment = useExperiment(experimentId);
  const corpus = useCorpusSamples();
  const [filter, setFilter] = useState("");

  const samples = useMemo(
    () => (corpus.data ?? []).filter((s) => s.experiment_id === experimentId),
    [corpus.data, experimentId],
  );
  const visible = useMemo(() => {
    const f = filter.trim().toLowerCase();
    return f ? samples.filter((s) => s.name.toLowerCase().includes(f)) : samples;
  }, [samples, filter]);

  const requested = sampleParam === undefined ? undefined : Number(sampleParam);
  const sample = requested === undefined ? samples[0] : samples.find((s) => s.id === requested);
  const sampleId = sample?.id;
  const goSample = useCallback(
    (sid: number) => navigate(`/experiments/${experimentId}/workspace/${sid}`, { replace: true }),
    [navigate, experimentId],
  );

  const exposures = useExposures(sampleId);
  const [exposureId, setExposureId] = useState<number | undefined>(linkedExposure);
  useEffect(() => setExposureId(linkedExposure), [sampleId, linkedExposure]);
  const exposureList: Exposure[] = useMemo(() => exposures.data ?? [], [exposures.data]);
  const exposure = exposureList.find((e) => e.id === exposureId) ?? representativeExposure(exposureList);

  const trace = useTrace(exposure?.id);
  const peaks = usePeaks(exposure?.id);
  const indices = useIndices(exposure?.id);
  const assignment = useAssignment(exposure?.id);
  const candidates = useMemo(() => rankCandidates(indices.data), [indices.data]);
  const assigned = useMemo(() => new Set(assignment.data?.members ?? []), [assignment.data]);
  const [hoverId, setHoverId] = useState<number | undefined>(undefined);

  const addPhase = useAddAssignmentPhase(exposure?.id ?? 0);
  const removePhase = useRemoveAssignmentPhase(exposure?.id ?? 0);
  const setCallState = useSetAssignmentState(exposure?.id ?? 0);
  const isFormFactor = assignment.data?.state === "form_factor";
  const selectExposure = useSelectExposure(sampleId ?? 0);
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
  const subRef = subChoice !== undefined && refs.includes(subChoice) ? subChoice : undefined;
  const [refsOn, setRefsOn] = useState(true);
  const [picking, setPicking] = useState(false);
  const [helpOpen, setHelpOpen] = useState(false);
  const [paletteOpen, setPaletteOpen] = useState(false);
  const theme = useWorkspaceTheme();
  const [draft, setDraft] = useState<CustomDraft | null>(null);
  const [customError, setCustomError] = useState<string | null>(null);
  useEffect(() => setDraft(null), [exposure?.id]);
  const draftOpen = draft !== null;
  useEffect(() => setCustomError(null), [draftOpen]);
  const commitCustom = useCommitCustomIndex(exposure?.id ?? 0);
  const compareTraces = useCompareTraces(refs);

  const onClickPeak = useCallback(
    (peakId: number) => {
      const action = peakClick(peaks.data ?? [], peakId);
      if (!action) return;
      if (action.kind === "remove") removePeak.mutate(peakId);
      else setExcluded.mutate({ peakId, excluded: action.excluded });
    },
    [peaks.data, removePeak, setExcluded],
  );

  const removeCompare = useCallback((rid: number) => {
    setRefs((r) => r.filter((x) => x !== rid));
    setSubRef((cur) => (cur === rid ? undefined : cur));
  }, []);

  const toggleCandidate = useCallback(
    (ix: IndexEntry) => {
      if (!exposure || exposure.status === "rejected") return;
      if (assigned.has(ix.id)) {
        removePhase.mutate(ix.id);
        return;
      }
      if (exposure.selected) {
        addPhase.mutate(ix.id);
        return;
      }
      selectExposure.mutate(exposure.id, { onSuccess: () => addPhase.mutate(ix.id) });
    },
    [exposure, assigned, addPhase, removePhase, selectExposure],
  );

  const toggleFormFactor = useCallback(() => {
    if (!exposure || exposure.status === "rejected") return;
    const next = isFormFactor ? "null" : "form_factor";
    if (exposure.selected) {
      setCallState.mutate(next);
      return;
    }
    selectExposure.mutate(exposure.id, { onSuccess: () => setCallState.mutate(next) });
  }, [exposure, isFormFactor, setCallState, selectExposure]);

  const toggleDrop = useCallback(
    (e: Exposure) => {
      setExposureId(e.id);
      setStatus.mutate({ exposureId: e.id, status: e.status === "rejected" ? null : "rejected" });
    },
    [setStatus],
  );

  const draftQs = useMemo(() => (draft ? draftPrediction(draft) : []), [draft]);
  const combs: CombLine[] = useMemo(() => {
    const shown = candidates.filter((c) => assigned.has(c.id) || c.id === hoverId);
    const lines = shown.map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase), emphasis: assigned.has(c.id) }));
    return draft ? [{ id: -1, qs: draftQs, color: phaseColor(draft.phase), emphasis: true }] : lines;
  }, [candidates, assigned, hoverId, draft, draftQs]);

  const sampleNames = useMemo(() => new Map(samples.map((s) => [s.id, s.name])), [samples]);
  const activeRefs = refs.filter((r) => r !== sampleId);
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
    return refs.flatMap((rid, i) => {
      const t = compareTraces.get(rid);
      const scale = scales.get(rid)?.scale;
      if (rid === sampleId || !t || scale === undefined) return [];
      const sc = scaledTrace(t, scale);
      return [{ id: rid, name: sampleNames.get(rid) ?? `Sample ${rid}`, q: sc.q, I: sc.I, dash: COMPARE_DASHES[i % COMPARE_DASHES.length]!, emphasis: rid === subRef }];
    });
  }, [mainTrace, refsOn, refs, compareTraces, scales, sampleId, sampleNames, subRef]);
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
  const pickable = useMemo(() => samples.filter((s) => s.id !== sampleId && !refs.includes(s.id)), [samples, sampleId, refs]);

  const activePeakQs = useMemo(() => (peaks.data ?? []).filter((p) => !p.excluded).map((p) => p.q), [peaks.data]);
  const canCall = Boolean(username && exposure && !exposureDropped);
  const draftBounds = draft ? latticeBounds(draft.phase, xExtent) : null;
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
    const commit = (): void => {
      commitCustom.mutate(c.phase, c.basis, c.ratios);
      setDraft(null);
    };
    if (exposure.selected) {
      commit();
      return;
    }
    setCustomError(null);
    selectExposure.mutate(exposure.id, {
      onSuccess: commit,
      onError: () => setCustomError("This frame could not be selected, so the index was not added. Try again."),
    });
  }, [draft, exposure, canCall, commitCustom, selectExposure]);
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
    (peakId: number) => {
      const p = (peaks.data ?? []).find((pk) => pk.id === peakId);
      if (p) snapDraft(p.q);
    },
    [peaks.data, snapDraft],
  );
  const draftFit = useMemo(() => (draft ? indexFit(draftQs, activePeakQs, xExtent) : null), [draft, draftQs, activePeakQs, xExtent]);

  const indexSeries: IndexSeries[] = useMemo(() => {
    if (!xExtent) return [];
    if (draft && draftFit) {
      return [{ key: "draft", phase: draft.phase, color: phaseColor(draft.phase), fit: draftFit, lattice: draft.lattice, r2: draftFit.r2 }];
    }
    const hovered = candidates.find((c) => c.id === hoverId);
    const shown = hovered ? [hovered] : candidates.filter((c) => assigned.has(c.id));
    return shown.map((c) => {
      const fit = storedFit(c, xExtent);
      return { key: String(c.id), phase: c.phase, color: phaseColor(c.phase), fit, lattice: c.lattice_d, r2: c.r_squared ?? fit.r2 };
    });
  }, [draft, draftFit, xExtent, candidates, hoverId, assigned]);

  useEffect(() => {
    const onKey = (ev: KeyboardEvent): void => {
      if ((ev.metaKey || ev.ctrlKey) && !ev.altKey && !ev.shiftKey && ev.key.toLowerCase() === "k") {
        ev.preventDefault();
        setHelpOpen(false);
        setPicking(false);
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
        case "j":
        case "J": {
          const nid = nextUnreviewed(visible, sampleId, ev.shiftKey ? -1 : 1);
          if (nid !== undefined) goSample(nid);
          break;
        }
        case "ArrowRight":
        case "ArrowLeft": {
          const eid = stepExposure(exposureList, exposure?.id, ev.key === "ArrowRight" ? 1 : -1);
          if (eid !== undefined) setExposureId(eid);
          break;
        }
        case "l":
        case "L":
          setView(nextView(view, ev.shiftKey ? -1 : 1));
          break;
        case "c":
        case "C":
          setPicking(true);
          break;
        case "i":
        case "I":
          openCustom();
          break;
        case "?":
          setHelpOpen(true);
          break;
        case "z":
        case "Z":
          setXDomain(null);
          break;
        case "r":
        case "R":
          setRefsOn((v) => !v);
          break;
        case "x":
        case "X":
          if (exposure && username) toggleDrop(exposure);
          break;
        case "0":
          if (username) toggleFormFactor();
          break;
        default: {
          const n = Number(ev.key);
          const c = Number.isInteger(n) && n >= 1 ? candidates[n - 1] : undefined;
          if (c && username) toggleCandidate(c);
          else handled = false;
        }
      }
      if (handled) ev.preventDefault();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [visible, sampleId, goSample, exposureList, exposure, toggleDrop, candidates, toggleCandidate, username, draft, addCustom, setDraftClamped, openCustom, toggleFormFactor, view, setView]);

  const reviewed = samples.filter((s) => sampleStatus(s).kind !== "unreviewed").length;
  const activePeaks = (peaks.data ?? []).filter((p) => !p.excluded).length;
  const exposureIndex = exposure ? exposureList.findIndex((e) => e.id === exposure.id) + 1 : 0;
  const callPhase = candidates.find((c) => assigned.has(c.id))?.phase ?? null;
  const emptyMessage = !Number.isInteger(experimentId)
    ? "This workspace link is not valid."
    : corpus.isLoading
      ? "Loading samples…"
      : samples.length === 0
        ? "No samples in this experiment."
        : "This sample is not in this experiment.";

  const commands: Command[] = [
    {
      id: "next-unreviewed",
      label: "Next unreviewed sample",
      shortcut: "J",
      run: () => {
        const nid = nextUnreviewed(visible, sampleId, 1);
        if (nid !== undefined) goSample(nid);
      },
    },
    { id: "compare", label: "Compare with a sample", shortcut: "C", run: () => setPicking(true) },
    { id: "refs", label: refsOn ? "Hide comparisons" : "Show comparisons", shortcut: "R", run: () => setRefsOn((v) => !v) },
    ...PLOT_VIEWS.map((o) => ({ id: `view-${o.view}`, label: `Plot ${o.name}${o.view === view ? ", current view" : ""}`, shortcut: "L", run: () => setView(o.view) })),
    { id: "zoom", label: "Reset zoom", shortcut: "Z", run: () => setXDomain(null) },
    { id: "theme", label: theme.isDark ? "Switch to light theme" : "Switch to dark theme", run: theme.toggle },
    ...(theme.theme === "system" ? [] : [{ id: "theme-system", label: "Use system theme", run: theme.followSystem }]),
    { id: "custom", label: "Custom index", shortcut: "I", run: openCustom },
    ...(username ? [{ id: "form-factor", label: isFormFactor ? "Clear no phase (form factor)" : "Mark no phase (form factor)", shortcut: "0", run: toggleFormFactor }] : []),
    { id: "help", label: "Keyboard shortcuts", shortcut: "?", run: () => setHelpOpen(true) },
    ...samples.map((s) => ({ id: `go-${s.id}`, label: `Go to ${s.name}`, hint: "Sample", run: () => goSample(s.id) })),
  ];

  return (
    <WorkspaceLayout
      theme={theme.theme}
      overlay={
        <>
          {helpOpen && <HelpOverlay sections={HELP} onClose={() => setHelpOpen(false)} />}
          {paletteOpen && <CommandPalette commands={commands} onClose={() => setPaletteOpen(false)} />}
        </>
      }
      header={
        <WorkspaceHeader
          experimentId={experimentId}
          experimentName={experiment.data?.name ?? `Experiment ${experimentId}`}
          reviewed={reviewed}
          total={samples.length}
          onHelp={() => setHelpOpen(true)}
        />
      }
      nav={<SampleList samples={visible} currentId={sampleId} filter={filter} onFilter={setFilter} onPick={goSample} />}
      main={
        sample ? (
          <>
            <Section title={sample.name} aside={exposure ? `Frame ${exposureIndex} of ${exposureList.length}` : undefined}>
              <div className="relative mb-2">
                <CompareBar
                  chips={activeRefs.map((rid) => ({ id: rid, name: sampleNames.get(rid) ?? `Sample ${rid}`, dash: COMPARE_DASHES[refs.indexOf(rid) % COMPARE_DASHES.length]!, scale: scales.get(rid) }))}
                  subtractId={subRef}
                  hidden={!refsOn}
                  view={view}
                  onView={setView}
                  onSubtract={(rid) => setSubRef((cur) => (cur === rid ? undefined : rid))}
                  onRemove={removeCompare}
                  onAdd={() => setPicking(true)}
                />
                {picking && (
                  <ComparePicker
                    samples={pickable}
                    onPick={(rid) => {
                      setRefs((r) => addCompare(r, rid));
                      if (subRef === undefined) setSubRef(rid);
                      setRefsOn(true);
                    }}
                    onClose={() => setPicking(false)}
                  />
                )}
              </div>
              <WorkspacePlot
                trace={trace.data}
                peaks={peaks.data ?? []}
                phase={callPhase}
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
            </Section>
            <Section title="Exposures">
              <ExposureStrip exposures={exposureList} currentId={exposure?.id} onPick={setExposureId} onToggleDrop={toggleDrop} />
            </Section>
          </>
        ) : (
          <Note>{emptyMessage}</Note>
        )
      }
      tools={
        <>
          {draft && draftBounds ? (
            <Section title="Custom index">
              <CustomIndexPanel
                draft={draft}
                min={draftBounds.min}
                max={draftBounds.max}
                landed={draftFit?.landed ?? 0}
                total={draftFit?.points.length ?? 0}
                disabled={!canCall || selectExposure.isPending}
                error={customError}
                onChange={setDraftClamped}
                onAdd={addCustom}
                onCancel={() => setDraft(null)}
              />
            </Section>
          ) : (
            <Section title="Phase">
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
          <Section title="Indexing">
            <IndexingPlot series={indexSeries} empty="Hover a candidate or tick one to see its fit." />
          </Section>
          {exposure && (
            <Section title="Exposure">
              <MetaList
                rows={[
                  ["File", exposure.filename ?? "—"],
                  ["Kind", exposure.kind === "file" ? "Raw" : exposure.kind === "averaged" ? "Averaged" : "Subtracted"],
                  ["Status", exposure.status === "rejected" ? "Dropped" : "Kept"],
                  ["Peaks", String(activePeaks)],
                ]}
              />
            </Section>
          )}
        </>
      }
    />
  );
}
