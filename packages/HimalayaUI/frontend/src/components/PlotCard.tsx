import { useCallback, useEffect, useMemo, useState } from "react";
import { useAppState } from "../state";
import {
  useExposures, useTrace, usePeaks, useIndices, useGroups,
  useAddPeak, useRemovePeak, useSetPeakExcluded,
  useExperiment, useSamples,
} from "../queries";
import { TraceViewer, interpolateI } from "./TraceViewer";
import { HintText } from "./ui";
import { phaseColor } from "../phases";
import type { IndexEntry, Peak, Trace } from "../api";

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
  const exposuresQ  = useExposures(activeSampleId, { excludeRejected: true });
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
  // Visible intensity-range (null = full data range). Auto-set by Fit features
  // and the per-exposure auto-fit; cleared by reset.
  const [yDomain, setYDomain] = useState<[number, number] | null>(null);
  // X-axis scale: log (SAXS convention) or linear.
  const [xType, setXType] = useState<"log" | "linear">("log");

  // Auto-pick first exposure when sample changes (or current choice is stale)
  useEffect(() => {
    const exposures = exposuresQ.data ?? [];
    if (exposures.length === 0) return;
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);

  // Reset the q-range when the sample or exposure changes — the previous
  // zoom almost never applies to a different trace.
  useEffect(() => { setXDomain(null); setYDomain(null); }, [activeExposureId]);

  // Compute a y-domain (and tightened x-domain when peaks exist) that focuses
  // on the diffraction features rather than the dominating beam decay.
  const computeFit = useCallback(
    (trace: Trace, peaks: Peak[]): {
      x: [number, number] | null;
      y: [number, number] | null;
    } => {
      if (trace.q.length < 2) return { x: null, y: null };

      if (peaks.length > 0) {
        const sortedQ = peaks.map((p) => p.q).sort((a, b) => a - b);
        const intensities = peaks
          .map((p) => interpolateI(p.q, trace))
          .filter((v) => Number.isFinite(v) && v > 0);
        if (intensities.length === 0) return { x: null, y: null };
        const lo = Math.min(...intensities);
        const hi = Math.max(...intensities);
        return {
          x: [sortedQ[0]! * 0.7, sortedQ[sortedQ.length - 1]! * 1.3],
          y: [lo / 5, hi * 5],
        };
      }

      // Fallback: drop the first 15% of x-range (the typical beam-decay
      // region) and use min/max intensity in the remainder.
      const startIdx = Math.floor(trace.q.length * 0.15);
      const tail = trace.I
        .slice(startIdx)
        .filter((v) => Number.isFinite(v) && v > 0);
      if (tail.length === 0) return { x: null, y: null };
      const lo = Math.min(...tail);
      const hi = Math.max(...tail);
      return { x: null, y: [lo / 2, hi * 2] };
    },
    [],
  );

  const fitFeatures = useCallback(() => {
    if (!traceQ.data) return;
    const fit = computeFit(traceQ.data, peaksQ.data ?? []);
    setXDomain(fit.x);
    setYDomain(fit.y);
  }, [traceQ.data, peaksQ.data, computeFit]);

  // Auto-fit once per exposure (re-fits when peaks finish loading too).
  useEffect(() => {
    if (!traceQ.data) return;
    const fit = computeFit(traceQ.data, peaksQ.data ?? []);
    setXDomain(fit.x);
    setYDomain(fit.y);
    // Intentionally only fires when the exposure changes or peaks stream in;
    // we don't want it to fight manual zoom on the same exposure.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeExposureId, peaksQ.data?.length, traceQ.data]);

  const resetDomain = useCallback(() => {
    setXDomain(null);
    setYDomain(null);
  }, []);

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
        yDomain={yDomain}
        xType={xType}
        onReset={resetDomain}
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
        xType={xType}
        onSetXType={setXType}
        onFitFeatures={fitFeatures}
        canFit={traceQ.data !== undefined}
      />
      <div className="relative flex-1 min-h-0">
        {body}
      </div>
      {activeExposureId !== undefined && traceQ.data && (
        <PlotLegend
          peaks={peaksQ.data ?? []}
          hoveredIndex={hoveredIndex}
        />
      )}
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
  xType: "log" | "linear";
  onSetXType: (t: "log" | "linear") => void;
  onFitFeatures: () => void;
  canFit: boolean;
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
  xType, onSetXType, onFitFeatures, canFit,
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
              <span className={hasExp ? "text-fg-muted" : "text-fg-muted italic"}>
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
      <div className="shrink-0 flex items-center gap-2">
        <button
          type="button"
          onClick={onFitFeatures}
          disabled={!canFit}
          data-testid="fit-features"
          title="Auto-zoom to peaks (or post-beam region)"
          className="text-[10.5px] px-1.5 py-0.5 rounded text-fg-dim hover:text-fg
                     hover:bg-bg-hover disabled:opacity-40 disabled:cursor-default
                     border border-transparent hover:border-border whitespace-nowrap"
        >
          fit features
        </button>
        <XScaleToggle xType={xType} onSetXType={onSetXType} />
        <QRange xDomain={xDomain} fullRange={fullRange} onXDomain={onXDomain} />
      </div>
    </div>
  );
}

interface XScaleToggleProps {
  xType: "log" | "linear";
  onSetXType: (t: "log" | "linear") => void;
}

function XScaleToggle({ xType, onSetXType }: XScaleToggleProps): JSX.Element {
  const btn = (val: "log" | "linear", label: string): JSX.Element => (
    <button
      type="button"
      onClick={() => onSetXType(val)}
      data-testid={`x-scale-${val}`}
      data-active={xType === val}
      className={[
        "text-[10.5px] px-1.5 py-0.5 transition-colors",
        xType === val
          ? "bg-bg-subtle text-fg"
          : "text-fg-dim hover:text-fg hover:bg-bg-hover",
      ].join(" ")}
    >
      {label}
    </button>
  );
  return (
    <span
      className="flex items-stretch border border-border rounded overflow-hidden"
      title="x-axis scale"
    >
      {btn("log", "log")}
      <span className="w-px bg-border" />
      {btn("linear", "lin")}
    </span>
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

// ── Plot legend ─────────────────────────────────────────────────────────────

interface PlotLegendProps {
  peaks: Peak[];
  hoveredIndex: IndexEntry | undefined;
}

function TriangleSvg({ color, opacity = 1 }: { color: string; opacity?: number }): JSX.Element {
  // Downward-pointing triangle matching TraceViewer geometry (hw=4, h=7)
  return (
    <svg width="10" height="8" viewBox="0 0 8 7" style={{ display: "block" }}>
      <polygon
        points="-4,0 4,0 0,7"
        transform="translate(4,0)"
        fill={color}
        fillOpacity={opacity}
        stroke={color}
        strokeOpacity={opacity}
        strokeWidth="0.5"
      />
    </svg>
  );
}

function TickLineSvg({ color }: { color: string }): JSX.Element {
  return (
    <svg width="6" height="12" viewBox="0 0 6 12" style={{ display: "block" }}>
      <line x1="3" y1="0" x2="3" y2="12" stroke={color} strokeWidth="1.5" strokeLinecap="round" />
    </svg>
  );
}

function LegendItem({
  symbol,
  label,
  style,
}: {
  symbol: JSX.Element;
  label: string;
  style?: React.CSSProperties;
}): JSX.Element {
  return (
    <span className="inline-flex items-center gap-1.5 whitespace-nowrap" style={style}>
      {symbol}
      {label}
    </span>
  );
}

function PlotLegend({ peaks, hoveredIndex }: PlotLegendProps): JSX.Element {
  const hasManualPeaks   = peaks.some((p) => p.source === "manual");
  const hasExcludedPeaks = peaks.some((p) => p.excluded);
  return (
    <div className="flex items-center gap-4 px-4 py-1.5 border-t border-border-soft
                    text-[10.5px] font-mono text-fg-dim flex-wrap">
      <LegendItem symbol={<TriangleSvg color="var(--color-accent)" />} label="auto peak" />
      {hasManualPeaks && (
        <LegendItem symbol={<TriangleSvg color="var(--color-peak-manual)" />} label="manual peak" />
      )}
      {hasExcludedPeaks && (
        <LegendItem symbol={<TriangleSvg color="var(--color-accent)" opacity={0.3} />} label="excluded" />
      )}
      {hoveredIndex && (
        <LegendItem
          symbol={<TickLineSvg color={phaseColor(hoveredIndex.phase)} />}
          label={`predicted ${hoveredIndex.phase}`}
          style={{ color: phaseColor(hoveredIndex.phase) }}
        />
      )}
    </div>
  );
}
