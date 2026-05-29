import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import {
  useTrace, usePeaks, useIndices, useGroups,
  useAddPeak, useRemovePeak, useSetPeakExcluded,
  useExperiment, useSamples, useExposures,
} from "../queries";
import { useAutoPickExposure, noUsableExposureState } from "../hooks/useAutoPickExposure";
import { TraceViewer } from "./TraceViewer";
import { NoUsableExposureNotice } from "./NoUsableExposureNotice";
import { Card, HintText } from "./ui";
import { SegmentedControl } from "./ui/SegmentedControl";
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
  const exposuresQ  = useExposures(activeSampleId);

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
  // R3-F03: "place a peak" armed mode. Presentational only — empty-area click
  // on the trace already adds a manual peak unconditionally; arming carries the
  // visible affordance + the `+ Peak` toggle state. Reset on exposure change.
  const [addArmed, setAddArmed] = useState(false);

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
  useEffect(() => { setXDomain(null); setYDomain(null); setAddArmed(false); }, [activeExposureId]);

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

  // R4 (L-11): when hovering a candidate that is NOT yet in the active set and
  // competes with an active phase (claims an overlapping peak), the peaks that
  // active phase would lose on the swap are dimmed in the trace so the cost is
  // visible before the click. Mirrors the mockup `setPreview` "losing" set.
  const losingPeakIds = useMemo(() => {
    if (!hoveredIndex) return undefined;
    const alreadyActive = activeGroupIndices.some((ix) => ix.id === hoveredIndex.id);
    if (alreadyActive) return undefined;
    const claim = new Set(hoveredIndex.peaks.map((p) => p.peak_id));
    const losing = new Set<number>();
    for (const active of activeGroupIndices) {
      const activePeakIds = active.peaks.map((p) => p.peak_id);
      const overlaps = activePeakIds.some((id) => claim.has(id));
      if (!overlaps) continue; // independent phase → coexists, nothing lost
      for (const id of activePeakIds) if (!claim.has(id)) losing.add(id);
    }
    return losing.size > 0 ? losing : undefined;
  }, [hoveredIndex, activeGroupIndices]);

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

  // Distinguish "no usable exposure" from a transient trace miss. When the
  // sample's exposures have loaded and none are acceptable (all rejected, or
  // none exist), activeExposureId stays undefined (useAutoPickExposure bails),
  // so the trace query never resolves. The M1 corpus→focus doors can land here
  // — surface a real empty state with a path back to the loupe rather than the
  // generic "No trace data available" copy, which reads as a load failure.
  const { noUsable: noUsableExposure, allRejected: allExposuresRejected } =
    noUsableExposureState(exposuresQ.data);

  const body = (() => {
    if (activeSampleId === undefined) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>Pick a sample to see its trace.</HintText>
        </div>
      );
    }
    if (noUsableExposure) {
      return (
        <NoUsableExposureNotice
          sampleId={activeSampleId}
          allRejected={allExposuresRejected}
        />
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
        losingPeakIds={losingPeakIds}
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
        addArmed={addArmed}
        onToggleAddArmed={() => setAddArmed((v) => !v)}
        {...(experimentQ.data?.q_units ? { qUnits: experimentQ.data.q_units } : {})}
      />
    );
  })();

  const titleStep: "experiment" | "sample" =
    activeExperimentId === undefined ? "experiment" : "sample";

  // R3-N3 (#209): the trace plate is the hero — the one elevated object in
  // The Print's "flat-except-the-plate" rule (DESIGN.md §Elevation). In the
  // focus variant it is the single elevated object (Card elevated). Non-focus
  // consumers stay an unstyled flow container (no plate, no bg) — gated on
  // headerSlot exactly as before.
  const layout = "flex flex-col h-full min-h-0 overflow-hidden";
  const inner = (
    <>
      <TitleStrip
        {...(headerSlot ? { headerSlot } : {})}
        experimentName={experimentName}
        sampleName={sampleName}
        onTitleClick={() => openNavModal(titleStep)}
        xType={xType}
        onSetXType={setXType}
        onFitFeatures={fitFeatures}
        canFit={traceQ.data !== undefined}
        isZoomed={xDomain !== null}
        addArmed={addArmed}
        onToggleAddPeak={() => setAddArmed((v) => !v)}
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
    </>
  );
  return headerSlot ? (
    <Card elevated data-testid="plot-card" className={layout}>
      {inner}
    </Card>
  ) : (
    <div data-testid="plot-card" className={layout}>
      {inner}
    </div>
  );
}

interface TitleStripProps {
  /** Focus-variant header (R3 / #226); replaces the picker when supplied. */
  headerSlot?:    JSX.Element;
  experimentName: string | undefined;
  sampleName:     string | undefined;
  onTitleClick:   () => void;
  xType: "log" | "linear";
  onSetXType: (t: "log" | "linear") => void;
  onFitFeatures: () => void;
  canFit: boolean;
  /** U-4: visible x-range is off the full trace (PlotCard's `xDomain !== null`). */
  isZoomed: boolean;
  /** R3-F03: `+ Peak` armed state + toggle. */
  addArmed: boolean;
  onToggleAddPeak: () => void;
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
  headerSlot, experimentName, sampleName, onTitleClick,
  xType, onSetXType, onFitFeatures, canFit, isZoomed, addArmed, onToggleAddPeak,
  exportSpec, exportFilenameStem, exportDisabled,
}: TitleStripProps): JSX.Element {
  const hasExp    = experimentName !== undefined;
  const hasSample = sampleName     !== undefined;
  // R3-N1 (#209): when `headerSlot` is supplied (the focus variant), use the
  // slotted card-header — drops the 56px clamp + bottom hairline that crushed
  // the focus header's 3-row stack flush against the edge. The legacy Index
  // path keeps the base `card-header` so other cards (PhasePanel, IndicesCard)
  // stay aligned with it.
  const headerClass = headerSlot
    ? "card-header card-header--slotted justify-between gap-3"
    : "card-header justify-between gap-3";
  return (
    <div
      data-testid="plot-stat-strip"
      data-variant={headerSlot ? "slotted" : "default"}
      className={headerClass}
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
                   hover:bg-paper-sunk transition-colors
                   focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent"
      >
        <span className="text-title tracking-tight truncate
                         max-w-[44ch]">
          {hasExp || hasSample ? (
            <>
              <span className={hasExp ? "text-ink-soft" : "text-ink-soft italic"}>
                {experimentName ?? "pick an experiment"}
              </span>
              <span className="text-ink-faint mx-1.5">·</span>
              <span className={hasSample ? "text-ink" : "text-ink-soft italic"}>
                {sampleName ?? "pick a sample"}
              </span>
            </>
          ) : (
            <span className="text-ink-soft italic">pick an experiment</span>
          )}
        </span>
        <span className="flex items-center gap-1.5 text-xs text-ink-faint leading-tight">
          <span>click to change</span>
          <span className="text-ink-faint/60">·</span>
          <kbd className="text-xs text-ink-faint
                          border border-hair-strong rounded px-1 leading-none py-px">/</kbd>
        </span>
      </button>
      )}
      <div className="shrink-0 flex items-center gap-2">
        {/* U-4: zoom-state indicator — only shows when the trace is off its
            full q-range; clicking it auto-fits (pairs with Auto-fit below). */}
        <ZoomIndicator zoomed={isZoomed} onReset={onFitFeatures} />
        <SegmentedControl<"log" | "linear">
          aria-label="x-axis scale"
          role="group"
          variant="bordered"
          options={[
            { value: "log", label: "log", testId: "x-scale-log" },
            { value: "linear", label: "lin", testId: "x-scale-linear" },
          ]}
          value={xType}
          onChange={onSetXType}
        />
        {/* R3-F03: the mockup's `.tools` cluster — Auto-fit + `+ Peak` ghost
            buttons replacing the numeric q-range input pair (anti-reference:
            "legacy scientific software" toolbar). Both are ghost `tool-btn`s
            matching the segmented scale toggle's text-xs scale. */}
        <button
          type="button"
          onClick={onFitFeatures}
          disabled={!canFit}
          data-testid="tool-autofit"
          title="Auto-zoom to peaks (or post-beam region)"
          className="rounded-md border border-hair-strong bg-plate px-2.5 py-1
                     text-xs font-semibold text-ink hover:bg-paper-sunk
                     disabled:opacity-40 disabled:cursor-default whitespace-nowrap"
        >
          Auto-fit
        </button>
        <button
          type="button"
          onClick={onToggleAddPeak}
          data-testid="tool-add-peak"
          data-armed={addArmed ? "true" : "false"}
          aria-pressed={addArmed}
          title="Click the trace to place a peak"
          className={[
            "rounded-md border px-2.5 py-1 text-xs font-semibold whitespace-nowrap transition-colors",
            addArmed
              ? "bg-print-accent border-print-accent text-paper"
              : "border-hair-strong bg-plate text-ink hover:bg-paper-sunk",
          ].join(" ")}
        >
          + Peak
        </button>
        {/* Thin divider before the export cluster. */}
        <span className="w-px h-4 bg-hair-strong" aria-hidden="true" />
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


/**
 * ZoomIndicator (U-4) — a quiet terracotta ghost button shown in the TitleStrip
 * only when the trace is zoomed off its full q-range. Clicking it auto-fits
 * (the parent passes `onFitFeatures`), so the chart "knows it's zoomed and
 * offers to reset" — pairing with the Auto-fit tool button.
 */
export function ZoomIndicator(
  { zoomed, onReset }: { zoomed: boolean; onReset: () => void },
): JSX.Element | null {
  if (!zoomed) return null;
  return (
    <button
      type="button"
      data-testid="zoom-indicator"
      onClick={onReset}
      title="Reset to full q-range"
      className="text-meta text-print-accent hover:underline whitespace-nowrap"
    >
      zoomed · reset
    </button>
  );
}

export interface QNumInputProps {
  value: number;
  onCommit: (v: number) => void;
  testId: string;
}

/**
 * QNumInput — a focus-gated numeric q input. R3-F03 retired the always-visible
 * `<QRange>` numeric pair from the trace tools cluster (it read as legacy
 * scientific-software toolbar soup); this component is retained, exported, and
 * unit-tested for a future scale-toggle popover that folds numeric q-range
 * editing back in off the segmented control. Not currently mounted.
 */
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
      className="w-[70px] bg-paper border border-hair-strong rounded px-1 py-0.5
                 text-ink text-xs tabular-nums text-right
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
    <div className="flex items-center gap-4 px-4 py-1.5 border-t border-hair
                    font-mono text-xs text-ink-faint flex-wrap">
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
