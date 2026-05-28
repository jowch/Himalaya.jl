import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import {
  useTrace, usePeaks, useIndices, useGroups,
  useAddPeak, useRemovePeak, useSetPeakExcluded,
  useExperiment, useSamples,
} from "../queries";
import { useAutoPickExposure } from "../hooks/useAutoPickExposure";
import { TraceViewer } from "./TraceViewer";
import { HintText } from "./ui";
import { FigureExportControls } from "./FigureExportControls";
import { buildTraceExportSpec } from "../lib/figure-export/adapters/traceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
import { phaseColor } from "../phases";
import type { IndexEntry, Peak, Trace } from "../api";

const PLOT_CARD_FIXTURE_DATA = {
  trace: {
    q:     [0.05,0.07,0.09,0.12,0.15,0.18,0.22,0.27,0.32,0.38,
            0.44,0.51,0.58,0.65,0.73,0.81,0.89,0.97,1.05,1.13],
    I:     [1200,980,820,650,520,420,350,290,240,310,
            280,190,150,120,95,80,68,60,54,50],
    sigma: Array<number>(20).fill(5),
  },
  peaks: [
    { id:1, exposure_id:0, q:0.18, intensity:420, prominence:180,
      sharpness:2.1, source:"auto" as const, excluded:false },
    { id:2, exposure_id:0, q:0.32, intensity:310, prominence:140,
      sharpness:1.8, source:"auto" as const, excluded:false },
    { id:3, exposure_id:0, q:0.51, intensity:190, prominence:100,
      sharpness:1.5, source:"auto" as const, excluded:false },
  ],
};

const PLOT_CARD_FIXTURE = (
  <TraceViewer
    trace={PLOT_CARD_FIXTURE_DATA.trace}
    peaks={PLOT_CARD_FIXTURE_DATA.peaks}
    activeGroupIndices={[]}
    hoveredIndex={undefined}
    onAddPeak={() => {}}
    onRemovePeak={() => {}}
    onTogglePeakExclusion={() => {}}
    xDomain={null}
    onXDomain={() => {}}
    yDomain={null}
    xType="log"
    onReset={() => {}}
  />
);

/**
 * PlotCard — center card on the Index page. Wraps the TraceViewer, a Miller-plot
 * inset anchored top-right, and a stat strip across the top.
 *
 * Auto-selects the first exposure when a sample is chosen (per the redesign plan
 * — exposure selection UI is deferred to a future triage page).
 */
export interface PlotCardProps {
  /**
   * Opt-in focus-variant header (R3 / #226). When supplied, the title strip
   * renders this node on the left in place of the legacy experiment-picker
   * button — dropping the `onTitleClick → openNavModal` affordance that, in
   * the focus context (route seeds only the sample), would otherwise fall to
   * the "pick an experiment" placeholder. The right-side q-controls / export
   * cluster is unchanged. Prop-less PlotCard keeps the picker (Index page).
   */
  headerSlot?: JSX.Element;
}

export function PlotCard({ headerSlot }: PlotCardProps = {}): JSX.Element {
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const activeSampleId     = useAppState((s) => s.activeSampleId);
  const activeExposureId   = useAppState((s) => s.activeExposureId);
  const hoveredIndexId     = useAppState((s) => s.hoveredIndexId);
  const hoveredPeakId      = useAppState((s) => s.hoveredPeakId);
  const hoveredQ           = useAppState((s) => s.hoveredQ);
  const setHoveredQ        = useAppState((s) => s.setHoveredQ);
  const openNavModal       = useAppState((s) => s.openNavModal);

  const experimentQ = useExperiment(activeExperimentId ?? 0);
  const samplesQ    = useSamples(activeExperimentId ?? 0);
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
    ? (sampleObj.display_name || sampleObj.name || `Sample ${activeSampleId}`)
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

  // Auto-pick the active exposure when the active sample changes — see the
  // hook for the full priority order. The regression test for issue #118
  // (test/sampleSwitchKeypress.test.tsx) subscribes to the same hook so
  // there's no manual mirror to drift.
  useAutoPickExposure(activeSampleId);

  // q-link (#180): stable adapter so the TraceViewer cursor effect (dep
  // [trace, peaks, onHoverQ]) doesn't tear down + re-add its mousemove
  // listeners every render. `setHoveredQ` is a stable Zustand action ref.
  const handleHoverQ = useCallback(
    (q: number | null) => setHoveredQ(q ?? undefined),
    [setHoveredQ],
  );

  // Reset the q-range when the sample or exposure changes — the previous
  // zoom almost never applies to a different trace.
  useEffect(() => { setXDomain(null); setYDomain(null); }, [activeExposureId]);

  // Compute a y-domain (and tightened x-domain when peaks exist) that focuses
  // on the diffraction features rather than the dominating beam decay.
  //
  // Y-floor strategy: derive from a low percentile of *positive* intensities
  // inside the visible x-window. This is the 1D analogue of the percentile
  // clip used on 2D images — it discards beamstop / dead-pixel zeros that
  // would otherwise drag the floor toward zero, while keeping the genuine
  // low-signal trace between peaks (so e.g. AgBe doesn't get clipped at the
  // bottom).
  const computeFit = useCallback(
    (trace: Trace, peaks: Peak[]): {
      x: [number, number] | null;
      y: [number, number] | null;
    } => {
      if (trace.q.length < 2) return { x: null, y: null };

      // 1. Pick the x-window. With peaks, bracket them; otherwise drop the
      //    first 15% of q (typical beam-decay region).
      let xMin = trace.q[0]!;
      let xMax = trace.q[trace.q.length - 1]!;
      let xResult: [number, number] | null = null;

      if (peaks.length > 0) {
        const sortedQ = peaks.map((p) => p.q).sort((a, b) => a - b);
        xMin = sortedQ[0]! * 0.7;
        xMax = sortedQ[sortedQ.length - 1]! * 1.3;
        xResult = [xMin, xMax];
      } else {
        const startIdx = Math.floor(trace.q.length * 0.15);
        xMin = trace.q[startIdx] ?? xMin;
      }

      // 2. Floor: 1st percentile of positive intensities WITHIN the focused
      //    x-window — close to the actual trace minimum so the low-signal
      //    tail is visible, but robust against single dead-pixel outliers.
      //    Filtering to the window means the floor reflects where the
      //    diffraction signal lives, not the beam region.
      const visibleI: number[] = [];
      for (let i = 0; i < trace.q.length; i++) {
        const q = trace.q[i]!;
        if (q < xMin || q > xMax) continue;
        const v = trace.I[i]!;
        if (Number.isFinite(v) && v > 0) visibleI.push(v);
      }
      if (visibleI.length === 0) return { x: xResult, y: null };
      visibleI.sort((a, b) => a - b);
      const p01 = visibleI[Math.floor(visibleI.length * 0.01)]!;

      // 3. Ceiling: actual max of the FULL trace (positive only). Keeping the
      //    full top end visible preserves relative-magnitude context — the
      //    user can see how much brighter the beam is than the peaks without
      //    having to reset the zoom. Only the floor is "cropped".
      let fullMax = 0;
      for (const v of trace.I) {
        if (Number.isFinite(v) && v > fullMax) fullMax = v;
      }
      if (fullMax <= 0) return { x: xResult, y: null };

      // Clamp the floor at fullMax / 1e5 so a stray near-zero pixel can't
      // blow the y-range past five log decades (which crushes the upper
      // signal into a tiny strip).
      const floor = Math.max(p01 * 0.5, fullMax / 1e5);
      return { x: xResult, y: [floor, fullMax * 1.2] };
    },
    [],
  );

  const fitFeatures = useCallback(() => {
    if (!traceQ.data) return;
    const fit = computeFit(traceQ.data, peaksQ.data ?? []);
    setXDomain(fit.x);
    setYDomain(fit.y);
  }, [traceQ.data, peaksQ.data, computeFit]);

  // Auto-fit once per exposure: trigger on the first render where both trace
  // and peaks are loaded for the active exposure. Subsequent peak mutations
  // (manual add/remove) don't re-fit — that would clobber the user's zoom.
  const fittedExposureRef = useRef<number | null>(null);
  useEffect(() => {
    if (activeExposureId === undefined) return;
    if (!traceQ.data || !peaksQ.data) return;
    if (fittedExposureRef.current === activeExposureId) return;
    fittedExposureRef.current = activeExposureId;
    const fit = computeFit(traceQ.data, peaksQ.data);
    setXDomain(fit.x);
    setYDomain(fit.y);
  }, [activeExposureId, traceQ.data, peaksQ.data, computeFit]);

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

  // Figure export (spec: 2026-05-08-figure-export-design.md).
  const exposureLabel = activeExposureId !== undefined
    ? `Exposure ${activeExposureId}`
    : "";
  const filenameStem = `himalaya-trace-${
    slugifyForFilename(experimentName ?? "")
  }-${
    slugifyForFilename(sampleName ?? "")
  }-${
    slugifyForFilename(exposureLabel)
  }`;

  const exportSpec = useCallback(() => {
    if (!traceQ.data || !peaksQ.data) {
      throw new Error("FigureExportControls: parent disabled-gate violated");
    }
    return buildTraceExportSpec({
      trace: traceQ.data,
      peaks: peaksQ.data,
      activeGroupIndices,
      experimentName: experimentName ?? "",
      sampleName: sampleName ?? "",
      exposureLabel,
      xDomain,
      yDomain,
      xType,
      ...(experimentQ.data?.q_units ? { qUnits: experimentQ.data.q_units } : {}),
    });
  }, [
    traceQ.data, peaksQ.data, activeGroupIndices,
    experimentName, sampleName, exposureLabel,
    xDomain, yDomain, xType, experimentQ.data?.q_units,
  ]);

  const exportDisabled = !traceQ.data || !peaksQ.data;

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
        hoveredPeakId={hoveredPeakId}
        hoveredQ={hoveredQ}
        onHoverQ={handleHoverQ}
        onAddPeak={(q) => addPeak.mutate(q)}
        onRemovePeak={(peakId) => removePeak.mutate(peakId)}
        onTogglePeakExclusion={(peakId, excluded) =>
          setPeakExcl.mutate({ peakId, excluded })}
        xDomain={xDomain}
        onXDomain={setXDomain}
        yDomain={yDomain}
        xType={xType}
        onReset={resetDomain}
        {...(experimentQ.data?.q_units ? { qUnits: experimentQ.data.q_units } : {})}
      />
    );
  })();

  const titleStep: "experiment" | "sample" =
    activeExperimentId === undefined ? "experiment" : "sample";

  return (
    <div data-testid="plot-card" className="flex flex-col h-full min-h-0 overflow-hidden">
      <TitleStrip
        headerSlot={headerSlot}
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
        exportSpec={exportSpec}
        exportFilenameStem={filenameStem}
        exportDisabled={exportDisabled}
      />
      <div className="relative flex-1 min-h-0">
        <Skeleton
          name="plot-card"
          className="h-full w-full"
          loading={activeExposureId !== undefined && (traceQ.isLoading || peaksQ.isLoading)}
          stagger={50}
          transition={200}
          fixture={PLOT_CARD_FIXTURE}
          fallback={<div className="flex-1 flex items-center justify-center"><HintText>Loading trace…</HintText></div>}
        >
          {body}
        </Skeleton>
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
  /** Focus-variant header (R3 / #226); replaces the picker when supplied. */
  headerSlot?:    JSX.Element;
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
  exportSpec: () => import("../lib/figure-export/types").ExportSpec;
  exportFilenameStem: string;
  exportDisabled: boolean;
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
  headerSlot, experimentName, sampleName, onTitleClick, xDomain, fullRange, onXDomain,
  xType, onSetXType, onFitFeatures, canFit,
  exportSpec, exportFilenameStem, exportDisabled,
}: TitleStripProps): JSX.Element {
  const hasExp    = experimentName !== undefined;
  const hasSample = sampleName     !== undefined;
  return (
    <div
      data-testid="plot-stat-strip"
      className="card-header justify-between gap-3"
    >
      {/* Focus variant (R3 / #226): render the supplied header and drop the
          experiment-picker button entirely. The picker only makes sense on the
          Index page, where both experiment and sample are picked globally. */}
      {headerSlot ?? (
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
        <span className="text-title tracking-tight truncate
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
        <span className="flex items-center gap-1.5 text-xs text-fg-dim leading-tight">
          <span>click to change</span>
          <span className="text-fg-dim/60">·</span>
          <kbd className="text-xs text-fg-dim
                          border border-border rounded px-1 leading-none py-px">/</kbd>
        </span>
      </button>
      )}
      <div className="shrink-0 flex items-center gap-2">
        <button
          type="button"
          onClick={onFitFeatures}
          disabled={!canFit}
          data-testid="fit-features"
          title="Auto-zoom to peaks (or post-beam region)"
          className="text-xs px-1.5 py-0.5 rounded text-fg-dim hover:text-fg
                     hover:bg-bg-hover disabled:opacity-40 disabled:cursor-default
                     border border-transparent hover:border-border whitespace-nowrap"
        >
          fit features
        </button>
        <XScaleToggle xType={xType} onSetXType={onSetXType} />
        <QRange xDomain={xDomain} fullRange={fullRange} onXDomain={onXDomain} />
        {/* Thin divider before the export cluster. */}
        <span className="w-px h-4 bg-border" aria-hidden="true" />
        <FigureExportControls
          spec={exportSpec}
          filenameStem={exportFilenameStem}
          ariaContext="trace plot"
          disabled={exportDisabled}
        />
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
        "text-xs px-1.5 py-0.5 transition-colors",
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
      <span className="text-fg-dim uppercase tracking-wider text-xs">q</span>
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
                 text-fg text-xs tabular-nums text-right
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
                    font-mono text-xs text-fg-dim flex-wrap">
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
