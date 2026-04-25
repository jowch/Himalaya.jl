import { useEffect, useMemo, useState } from "react";
import { useAppState } from "../state";
import {
  useExposures, useTrace, usePeaks, useIndices, useGroups,
  useAddPeak, useRemovePeak, useSetPeakExcluded,
  useExperiment, useSamples,
} from "../queries";
import { TraceViewer } from "./TraceViewer";
import { HintText } from "./ui";

/**
 * PlotCard — center card on the Index page. Wraps the TraceViewer, a Miller-plot
 * inset anchored top-right, and a stat strip across the top.
 *
 * Auto-selects the first exposure when a sample is chosen (per the redesign plan
 * — exposure selection UI is deferred to a future triage page).
 */
export function PlotCard(): JSX.Element {
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const activeSampleId     = useAppState((s) => s.activeSampleId);
  const activeExposureId   = useAppState((s) => s.activeExposureId);
  const hoveredIndexId     = useAppState((s) => s.hoveredIndexId);
  const setActiveExposure  = useAppState((s) => s.setActiveExposure);
  const openNavModal       = useAppState((s) => s.openNavModal);

  const experimentQ = useExperiment(activeExperimentId ?? 0);
  const samplesQ    = useSamples(activeExperimentId ?? 0);
  const exposuresQ  = useExposures(activeSampleId);
  const traceQ      = useTrace(activeExposureId);
  const peaksQ      = usePeaks(activeExposureId);
  const indicesQ    = useIndices(activeExposureId);
  const groupsQ     = useGroups(activeExposureId);

  const experimentName = activeExperimentId !== undefined
    ? (experimentQ.data?.name ?? `Experiment ${activeExperimentId}`)
    : undefined;
  const sampleObj = activeSampleId !== undefined
    ? samplesQ.data?.find((s) => s.id === activeSampleId)
    : undefined;
  const sampleName = sampleObj
    ? (sampleObj.name ?? sampleObj.label ?? `Sample ${activeSampleId}`)
    : undefined;

  const addPeak       = useAddPeak(activeExposureId ?? 0);
  const removePeak    = useRemovePeak(activeExposureId ?? 0);
  const setPeakExcl   = useSetPeakExcluded(activeExposureId ?? 0);

  // Visible q-range (null = full trace). Shared between TraceViewer wheel-zoom
  // and the numeric inputs in StatStrip.
  const [xDomain, setXDomain] = useState<[number, number] | null>(null);

  // Auto-pick first exposure when sample changes (or current choice is stale)
  useEffect(() => {
    const exposures = exposuresQ.data ?? [];
    if (exposures.length === 0) return;
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);

  // Reset the q-range when the sample or exposure changes — the previous
  // zoom almost never applies to a different trace.
  useEffect(() => { setXDomain(null); }, [activeExposureId]);

  const indices = indicesQ.data ?? [];
  const activeGroup = (groupsQ.data ?? []).find((g) => g.active);
  const activeGroupIndices = useMemo(
    () => (activeGroup?.members ?? [])
      .map((id) => indices.find((i) => i.id === id))
      .filter((i): i is NonNullable<typeof i> => i != null),
    [activeGroup, indices],
  );
  const hoveredIndex = hoveredIndexId != null
    ? indices.find((i) => i.id === hoveredIndexId)
    : undefined;

  const fullQRange: [number, number] | null = traceQ.data && traceQ.data.q.length > 0
    ? [traceQ.data.q[0]!, traceQ.data.q[traceQ.data.q.length - 1]!]
    : null;
  const effectiveDomain = xDomain ?? fullQRange;

  const body = (() => {
    if (activeSampleId === undefined) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>Pick a sample to see its trace.</HintText>
        </div>
      );
    }
    if (activeExposureId === undefined || traceQ.isPending || peaksQ.isPending) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>Loading trace…</HintText>
        </div>
      );
    }
    if (!traceQ.data || !peaksQ.data) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>No trace data available.</HintText>
        </div>
      );
    }
    return (
      <TraceViewer
        trace={traceQ.data}
        peaks={peaksQ.data}
        activeGroupIndices={activeGroupIndices}
        hoveredIndex={hoveredIndex}
        onAddPeak={(q) => addPeak.mutate(q)}
        onRemovePeak={(peakId) => removePeak.mutate(peakId)}
        onTogglePeakExclusion={(peakId, excluded) =>
          setPeakExcl.mutate({ peakId, excluded })}
        xDomain={xDomain}
        onXDomain={setXDomain}
      />
    );
  })();

  const titleStep: "experiment" | "sample" =
    activeExperimentId === undefined ? "experiment" : "sample";

  return (
    <div data-testid="plot-card" className="flex flex-col h-full min-h-0 overflow-hidden">
      <TitleStrip
        experimentName={experimentName}
        sampleName={sampleName}
        onTitleClick={() => openNavModal(titleStep)}
        xDomain={effectiveDomain}
        fullRange={fullQRange}
        onXDomain={setXDomain}
      />
      <div className="relative flex-1 min-h-0">
        {body}
      </div>
    </div>
  );
}

interface TitleStripProps {
  experimentName: string | undefined;
  sampleName:     string | undefined;
  onTitleClick:   () => void;
  xDomain:  [number, number] | null;
  fullRange: [number, number] | null;
  onXDomain: (d: [number, number] | null) => void;
}

/**
 * TitleStrip — top of the plot card.
 *
 * Layout:
 *   [ Experiment · Sample          ] [        q-range controls ]
 *   [ click to change · /          ]
 *
 * Two text rows so the strip is the same height as the PhasePanel header in
 * the right card — the two cards then read as a single horizontal band at
 * the top of the workspace.
 */
function TitleStrip({
  experimentName, sampleName, onTitleClick, xDomain, fullRange, onXDomain,
}: TitleStripProps): JSX.Element {
  const hasExp    = experimentName !== undefined;
  const hasSample = sampleName     !== undefined;
  return (
    <div
      data-testid="plot-stat-strip"
      className="card-header justify-between gap-3"
    >
      <button
        type="button"
        data-testid="plot-title"
        onClick={onTitleClick}
        title="Change experiment / sample"
        className="min-w-0 flex flex-col items-start text-left
                   -mx-1 px-1 py-0.5 rounded
                   hover:bg-bg-hover transition-colors
                   focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent"
      >
        <span className="text-[13px] font-semibold tracking-tight truncate
                         max-w-[44ch]">
          {hasExp || hasSample ? (
            <>
              <span className={hasExp ? "text-fg" : "text-fg-muted italic"}>
                {experimentName ?? "pick an experiment"}
              </span>
              <span className="text-fg-dim mx-1.5">·</span>
              <span className={hasSample ? "text-fg" : "text-fg-muted italic"}>
                {sampleName ?? "pick a sample"}
              </span>
            </>
          ) : (
            <span className="text-fg-muted italic">pick an experiment</span>
          )}
        </span>
        <span className="flex items-center gap-1.5 text-[11px] text-fg-dim leading-tight">
          <span>click to change</span>
          <span className="text-fg-dim/60">·</span>
          <kbd className="text-[10px] text-fg-dim
                          border border-border rounded px-1 leading-none py-px">/</kbd>
        </span>
      </button>
      <div className="shrink-0">
        <QRange xDomain={xDomain} fullRange={fullRange} onXDomain={onXDomain} />
      </div>
    </div>
  );
}

interface QRangeProps {
  xDomain: [number, number] | null;
  fullRange: [number, number] | null;
  onXDomain: (d: [number, number] | null) => void;
}

function QRange({ xDomain, fullRange, onXDomain }: QRangeProps): JSX.Element | null {
  if (!fullRange) return null;
  const [qmin, qmax] = xDomain ?? fullRange;
  const isFull = !xDomain;

  const commit = (nextMin: number, nextMax: number): void => {
    const lo = Math.max(fullRange[0], Math.min(nextMin, nextMax));
    const hi = Math.min(fullRange[1], Math.max(nextMin, nextMax));
    if (hi - lo < (fullRange[1] - fullRange[0]) * 1e-4) return;
    onXDomain([lo, hi]);
  };

  return (
    <span
      className="flex items-center gap-1.5 whitespace-nowrap"
      data-testid="q-range-controls"
    >
      <span className="text-fg-dim uppercase tracking-wider text-[9.5px]">q</span>
      <QNumInput
        value={qmin}
        onCommit={(v) => commit(v, qmax)}
        testId="q-range-min"
      />
      <span className="text-fg-dim">–</span>
      <QNumInput
        value={qmax}
        onCommit={(v) => commit(qmin, v)}
        testId="q-range-max"
      />
      <button
        type="button"
        onClick={() => onXDomain(null)}
        disabled={isFull}
        data-testid="q-range-reset"
        title="Reset q-range (double-click plot)"
        className="ml-1 px-1.5 py-0.5 rounded text-fg-dim hover:text-fg
                   hover:bg-bg-hover disabled:opacity-40 disabled:cursor-default
                   border border-transparent hover:border-border"
      >
        reset
      </button>
    </span>
  );
}

export interface QNumInputProps {
  value: number;
  onCommit: (v: number) => void;
  testId: string;
}

export function QNumInput({ value, onCommit, testId }: QNumInputProps): JSX.Element {
  const [draft, setDraft] = useState(value.toFixed(3));
  const [focused, setFocused] = useState(false);

  // Sync external value changes into the draft only when not actively editing.
  useEffect(() => {
    if (!focused) setDraft(value.toFixed(3));
  }, [value, focused]);

  return (
    <input
      type="number"
      step="0.001"
      value={draft}
      data-testid={testId}
      onChange={(e) => setDraft(e.currentTarget.value)}
      onFocus={() => setFocused(true)}
      onBlur={(e) => {
        setFocused(false);
        const n = parseFloat(e.currentTarget.value);
        if (Number.isFinite(n)) onCommit(n);
      }}
      onKeyDown={(e) => {
        if (e.key === "Enter") {
          const n = parseFloat((e.currentTarget as HTMLInputElement).value);
          if (Number.isFinite(n)) onCommit(n);
          (e.currentTarget as HTMLInputElement).blur();
        }
      }}
      className="w-[70px] bg-bg border border-border rounded px-1 py-0.5
                 text-fg text-[10.5px] tabular-nums text-right
                 outline-0 focus:border-accent"
    />
  );
}
