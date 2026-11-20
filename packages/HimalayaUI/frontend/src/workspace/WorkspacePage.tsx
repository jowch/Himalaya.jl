import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import type { Exposure, IndexEntry } from "../api";
import {
  useAddAssignmentPhase,
  useAssignment,
  useCorpusSamples,
  useExperiment,
  useExposures,
  useIndices,
  usePeaks,
  useRemoveAssignmentPhase,
  useSelectExposure,
  useSetExposureStatus,
  useTrace,
} from "../queries";
import { representativeExposure } from "../lib/exposure";
import { phaseColor } from "../phases";
import { useAppState } from "../state";
import { nextUnreviewed, rankCandidates, sampleStatus, stepExposure } from "./model";
import { WorkspaceHeader } from "./ui/WorkspaceHeader";
import { Note } from "./ui/Note";
import { WorkspaceLayout } from "./ui/WorkspaceLayout";
import { SampleList } from "./ui/SampleList";
import { WorkspacePlot, type CombLine } from "./ui/WorkspacePlot";
import { ExposureStrip } from "./ui/ExposureStrip";
import { CandidateList } from "./ui/CandidateList";
import { MetaList } from "./ui/MetaList";
import { Section } from "./ui/Section";

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
  const [exposureId, setExposureId] = useState<number | undefined>(undefined);
  useEffect(() => setExposureId(undefined), [sampleId]);
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
  const selectExposure = useSelectExposure(sampleId ?? 0);
  const setStatus = useSetExposureStatus(sampleId ?? 0);

  const exposureDropped = exposure?.status === "rejected";

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

  const toggleDrop = useCallback(
    (e: Exposure) => {
      setExposureId(e.id);
      setStatus.mutate({ exposureId: e.id, status: e.status === "rejected" ? null : "rejected" });
    },
    [setStatus],
  );

  const combs: CombLine[] = useMemo(() => {
    const shown = candidates.filter((c) => assigned.has(c.id) || c.id === hoverId);
    return shown.map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase), emphasis: assigned.has(c.id) }));
  }, [candidates, assigned, hoverId]);

  useEffect(() => {
    const onKey = (ev: KeyboardEvent): void => {
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
        case "x":
        case "X":
          if (exposure && username) toggleDrop(exposure);
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
  }, [visible, sampleId, goSample, exposureList, exposure, toggleDrop, candidates, toggleCandidate, username]);

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

  return (
    <WorkspaceLayout
      header={
        <WorkspaceHeader
          experimentId={experimentId}
          experimentName={experiment.data?.name ?? `Experiment ${experimentId}`}
          reviewed={reviewed}
          total={samples.length}
        />
      }
      nav={<SampleList samples={visible} currentId={sampleId} filter={filter} onFilter={setFilter} onPick={goSample} />}
      main={
        sample ? (
          <>
            <Section title={sample.name} aside={exposure ? `Frame ${exposureIndex} of ${exposureList.length}` : undefined}>
              <WorkspacePlot trace={trace.data} peaks={peaks.data ?? []} phase={callPhase} combs={combs} />
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
          <Section title="Phase">
            {exposureDropped && <Note>This frame is dropped. Restore it to call a phase.</Note>}
            <CandidateList
              candidates={candidates}
              assigned={assigned}
              disabled={exposureDropped || !username}
              onToggle={toggleCandidate}
              onHover={setHoverId}
            />
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
