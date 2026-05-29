import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useSeries, useMemberTraces, useMemberTracesLoading,
  useMemberExposures, useMemberSamples,
} from "../queries";
import { MultiTracePlot, offsetToBandFraction } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { GroupingModeToggle } from "../components/GroupingModeToggle";
import { AnnotationToggles } from "../components/AnnotationToggles";
import { ActiveBandProvider } from "../components/ActiveBandContext";
import { FigureExportControls } from "../components/FigureExportControls";
import { SeriesBuilderRail } from "../components/SeriesBuilderRail";
import { SeriesRecipeEditor } from "../components/SeriesRecipeEditor";
import { OffsetDock } from "../components/OffsetDock";
import { HintText } from "../components/ui";
import type { Representation } from "../components/RepresentationToggle";
import type { ScaleMode } from "../components/ScaleToggle";
import { resolveDisplayLabels } from "../lib/comparison/labels";
import { buildMultiTraceExportSpec } from "../lib/figure-export/adapters/multiTraceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
import { useAppState } from "../state";
import type { GroupingMode } from "../lib/comparison/coloring";
import type { Series, SeriesMember } from "../api";

/** Static skeleton for boneyard's headless capture: the plate area + a rail. */
const BUILDER_FIXTURE = (
  <div className="grid grid-cols-[1fr_336px] gap-0">
    <div className="m-4 rounded border border-hair bg-plate" style={{ aspectRatio: "10 / 3" }} />
    <div className="border-l border-hair p-4">
      <div className="h-4 w-1/2 rounded bg-paper-sunk" />
    </div>
  </div>
);

const GROUPING_MODES: readonly GroupingMode[] = ["bySample", "byPhase", "distinct"];

/**
 * Validate the series' persisted (intentionally-loose `string | null`)
 * view_grouping_mode against the GroupingMode union, falling back to
 * "bySample" (effectiveGroupingMode's hard default). A membership guard, not
 * a blind cast — a stored value outside the union can't slip through.
 */
function coerceGroupingMode(value: string | null): GroupingMode {
  return GROUPING_MODES.includes(value as GroupingMode)
    ? (value as GroupingMode)
    : "bySample";
}

/**
 * SeriesBuilderPage — the series builder visual surface at /series/:id
 * (#175 / I3.5a). Read-only: reads one series via useSeries(id) and composes
 * the existing MultiTracePlot render core. Mutations (recipe edits, plate
 * commit, permalink) are I3.5b — NOT here. Mounted under the CorpusShell
 * layout route, the destination the I3.3 folio card links to.
 *
 * URL-owned: the series id comes from the route param, never Zustand.
 */
export function SeriesBuilderPage(): JSX.Element {
  const { id: idParam } = useParams<{ id: string }>();
  const seriesId = Number(idParam);
  const query = useSeries(Number.isFinite(seriesId) ? seriesId : undefined);

  if (query.isError) {
    return (
      <div data-testid="series-builder-error" className="p-8 text-sm text-ink-soft">
        Could not load this series. It may have been deleted.
      </div>
    );
  }

  const s = query.data;

  // I3.5b — edit-mode gate. A draft for THIS series being active flips the
  // builder into edit mode (the recipe editor injects into the rail). The
  // Edit button seeds the draft; entering edit mode is draft-driven on the
  // bare /series/:id URL (no /edit segment), mirroring Compare.
  const seriesDraft = useAppState((st) => st.seriesDraft);
  const startSeriesDraft = useAppState((st) => st.startSeriesDraftFromSeries);
  const editing = s !== undefined && seriesDraft !== null && seriesDraft.id === s.id;

  return (
    <div data-testid="series-builder-page" className="flex h-full min-h-0 flex-col">
      <Skeleton
        name="series-builder"
        className="flex-1 min-h-0 flex flex-col"
        loading={query.isLoading}
        fixture={BUILDER_FIXTURE}
        fallback={
          <div className="flex-1 flex items-center justify-center">
            <HintText>Loading series…</HintText>
          </div>
        }
      >
        {/*
          R8-N1 (round-2): no outer page <header> here. The figure-plate
          kicker+title inside SeriesBuilderBody is the only header on this
          surface — see `series-builder.html:386-396`. The Edit button +
          "Editing" badge are re-homed into the plate kicker row so the
          "figure-as-plate" metaphor reads cleanly.
        */}
        {s && (
          <SeriesBuilderBody
            series={s}
            editing={editing}
            onStartEdit={() => startSeriesDraft(s)}
          />
        )}
      </Skeleton>
    </div>
  );
}

/**
 * Loaded body — composes the render core once `series` is present. Mirrors
 * Compare.tsx's review body wiring (the single source of truth for prop
 * shapes); the only differences are the Series input type and the local-state
 * groupingMode/xDomain (no draft, no Zustand comparison-keyed domain — a
 * read surface's coloring + pan/zoom are local UI concerns).
 */
function SeriesBuilderBody(
  { series: s, editing, onStartEdit }: {
    series: Series;
    editing: boolean;
    /** Seeds a fresh draft for this series — wired by the Edit button in the
     *  plate kicker row (R8-N1). Also wired to the autogroup card's "Adjust"
     *  affordance so both entry points share one action. */
    onStartEdit: () => void;
  },
): JSX.Element {
  // Members arrive sorted by display_order from the route; keep that order.
  const members: SeriesMember[] = s.members;
  const exposureIds = useMemo(
    () => members.flatMap((m) => (m.exposure_id !== null ? [m.exposure_id] : [])),
    [members],
  );

  // Hydrate per-member trace / exposure / sample rows via the existing hooks.
  // Each hook is called exactly once at top level (Rules of Hooks).
  const traces = useMemberTraces(exposureIds);
  const tracesLoading = useMemberTracesLoading(exposureIds);
  const exposures = useMemberExposures(exposureIds);
  const sampleIds = useMemo(() => {
    const ids = new Set<number>();
    for (const e of exposures.values()) ids.add(e.sample_id);
    return Array.from(ids).sort((a, b) => a - b);
  }, [exposures]);
  const samples = useMemberSamples(sampleIds);

  // Coloring mode: seed from the series' persisted view default, else
  // "bySample" — matching effectiveGroupingMode's hard default so the
  // builder's default coloring is identical to Compare's. The helper itself
  // is not reused: it is Comparison-typed and draft-aware, and this read
  // surface has no draft (persistence is a recipe edit → I3.5b). Local state.
  const [groupingMode, setGroupingMode] = useState<GroupingMode>(
    coerceGroupingMode(s.view_grouping_mode),
  );

  // Pan/zoom q-domain. MultiTracePlot requires non-optional xDomain +
  // onXDomain. The builder is a read surface, so its pan/zoom is local —
  // held in local useState rather than Zustand. null = full data range.
  // (Compare's per-comparison-id `compareXDomains` Zustand record was removed
  // with the Compare page in I5.3 #184; a series-id reuse would have
  // namespace-collided anyway.)
  const [xDomain, setXDomain] = useState<[number, number] | null>(null);

  // sampleIdFor: exact signature from MultiTracePlotProps / Compare.tsx —
  // `(m: SeriesMember) => number | null` (null, not a -1 sentinel).
  const sampleIdFor = useCallback(
    (m: SeriesMember): number | null => {
      if (m.exposure_id === null) return null;
      return exposures.get(m.exposure_id)?.sample_id ?? null;
    },
    [exposures],
  );

  // arg order: resolveDisplayLabels(members, exposures, samples).
  const displayLabelByMemberId = useMemo(
    () => resolveDisplayLabels(members, exposures, samples),
    [members, exposures, samples],
  );

  // Local UI state for the rail: collapse (full-bleed) + representation
  // (waterfall + heatmap, both live in the shared render core, #208).
  const [collapsed, setCollapsed] = useState(false);
  const [representation, setRepresentation] = useState<Representation>("waterfall");

  // Cross-trace peak-tracking layer (#208). Off by default; mockup
  // `series-builder.html:474-476` describes it as a Display option that
  // "links each Miller order across the traces that carry it — reads best
  // when peaks align." Local UI state; persisting it is a future view-pref.
  const [trackOn, setTrackOn] = useState(false);

  // Compose controls (R8 / B-F): the trace-offset slider and the log/linear
  // q-axis scale. Local UI state — a read surface's composition is a local
  // concern (persisting it is a recipe edit, deferred). Offset maps to the
  // plot's working-band fraction; scaleMode maps to xType.
  const [offset, setOffset] = useState(1.2);
  const [scaleMode, setScaleMode] = useState<ScaleMode>("log");
  const workingBandFraction = offsetToBandFraction(offset);

  // "Adjust" on the autogroup card enters edit mode — same affordance as the
  // plate kicker's Edit button. Both call `onStartEdit` so the seed action
  // lives in exactly one place.

  // Annotation toggles live in Zustand (shared with AnnotationToggles); the
  // export spec must reflect their current value so the figure matches the
  // on-screen plot.
  const showPeakTicks = useAppState((st) => st.showPeakTicks);
  const showPeakLabels = useAppState((st) => st.showPeakLabels);

  // Figure export — mirror Compare's spec thunk (evaluated at click time so
  // it captures fresh xDomain / toggles). No experiment scope on this surface,
  // so experimentName is omitted and the filename stem uses the series title.
  // representation + showCrossTraceTracking propagate so the exported PNG/SVG
  // matches the on-screen plot when heatmap or tracking is active (#251 r1).
  const exportFilenameStem = `himalaya-series-${slugifyForFilename(s.title || String(s.id))}`;
  const exportSpec = useCallback(() => buildMultiTraceExportSpec({
    members,
    traces,
    comparisonTitle: s.title,
    xDomain,
    showPeakTicks,
    showPeakLabels,
    groupingMode,
    sampleIdFor,
    displayLabelByMemberId,
    representation,
    showCrossTraceTracking: trackOn,
  }), [
    members, traces, s.title, xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
    representation, trackOn,
  ]);
  const exportDisabled =
    tracesLoading
    || members.length === 0
    || traces.size === 0
    || members.every((m) => m.exposure_id === null);

  // Track the plot column height so the gutter rows align with the y-bands
  // (both consumers share computeYBands). Mirrors Compare's ResizeObserver.
  const plotColRef = useRef<HTMLDivElement>(null);
  const [panelHeight, setPanelHeight] = useState(0);
  useEffect(() => {
    const el = plotColRef.current;
    if (!el) return;
    setPanelHeight(el.clientHeight);
    const obs = new ResizeObserver(() => {
      if (plotColRef.current) setPanelHeight(plotColRef.current.clientHeight);
    });
    obs.observe(el);
    return () => obs.disconnect();
  }, [tracesLoading]);

  return (
    <ActiveBandProvider>
      <div className="flex-1 min-h-0 flex flex-row" data-testid="series-builder-body">
        <div
          className="flex-1 min-h-0 overflow-auto flex flex-col items-center px-8 py-6"
          data-testid="series-builder-plot"
        >
          {/*
            Figure-as-plate (R8 / B-J): the figure is the printed plate —
            white plate, hair border, soft shadow, centered, widening to
            1336px when the rail collapses to full-bleed.

            After R8-N1 (round-2) the plate is the SINGLE header on this
            surface (mockup `series-builder.html:386-396`). The plate renders
            in both empty and populated states so the title + Edit
            affordance stay reachable — adding the first sample is the
            recipe-edit flow that populates the plot.
          */}
          <div
            data-testid="series-builder-plate"
            className={`w-full ${collapsed ? "max-w-[1336px]" : "max-w-[1180px]"} rounded border border-hair bg-plate p-8 shadow-[0_1px_0_rgba(255,255,255,.6)_inset,0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-20px_rgba(60,52,40,.22)] transition-[max-width] duration-200`}
          >
              {/*
                Kicker tag-row (R8 / B-H): terracotta "Series" + facet tags.
                After R8-N1 (round-2) the Edit button / "Editing" badge sit
                inline at the row's right edge — the page no longer ships an
                outer header strip, so this is where the affordance lives.
                Mockup `series-builder.html:386-396` shows only the plate
                header; the editing affordance is page-only chrome and lives
                with the kicker as the lightest re-home.
              */}
              <div
                className="mb-2 flex items-baseline justify-between gap-3"
                data-testid="fig-tags"
              >
                <div className="flex items-baseline gap-3">
                  <span className="text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">
                    Series
                  </span>
                  <div className="flex gap-1.5">
                    <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                      {members.length} {members.length === 1 ? "sample" : "samples"}
                    </span>
                    <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                      {scaleMode === "log" ? "log q" : "linear q"}
                    </span>
                    <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                      {representation === "heatmap" ? "intensity map" : "offset waterfall"}
                    </span>
                    {/* tracked tag — mockup `series-builder.html:393` (#251 r1 / N3) */}
                    {trackOn && (
                      <span
                        data-testid="fig-tag-track"
                        className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint"
                      >
                        tracked
                      </span>
                    )}
                  </div>
                </div>
                {editing ? (
                  <span data-testid="series-builder-editing-badge" className="text-[10.5px] uppercase tracking-wide text-ink-faint">
                    Editing
                  </span>
                ) : (
                  <button
                    type="button"
                    data-testid="series-builder-edit"
                    onClick={onStartEdit}
                    className="rounded border border-hair px-2 py-0.5 text-[11px] text-ink hover:bg-paper-sunk"
                  >
                    Edit
                  </button>
                )}
              </div>
              <h1 className="text-display font-medium text-ink">{s.title || "Untitled series"}</h1>
              {s.description && (
                <p className="mt-2 max-w-[64ch] text-sm text-ink-soft" data-testid="fig-sub">
                  {s.description}
                </p>
              )}

              {members.length === 0 ? (
                // Empty-state placeholder lives INSIDE the plate so the kicker
                // + title + Edit affordance remain visible. Adding the first
                // sample is exactly the I3.5b recipe-edit flow; gating the
                // editor out of the empty case would lock that flow out.
                <div
                  data-testid="series-builder-empty"
                  className="mt-6 grid place-items-center gap-3 rounded border border-dashed border-hair px-6 py-16 text-sm text-ink-faint"
                >
                  <p>This series has no members yet.</p>
                  {!editing && (
                    <button
                      type="button"
                      data-testid="series-builder-empty-cta"
                      onClick={onStartEdit}
                      className="font-medium text-print-accent hover:underline focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
                    >
                      Add the first sample
                    </button>
                  )}
                </div>
              ) : (
                <>
                  <div className="mt-4 flex items-center gap-3" data-testid="series-builder-controls">
                    <GroupingModeToggle value={groupingMode} onChange={setGroupingMode} />
                    <AnnotationToggles />
                  </div>
                  <div className="mt-3 flex min-h-0 flex-row gap-2" style={{ height: "60vh" }}>
                    <div ref={plotColRef} className="relative min-w-0 flex-1">
                      {representation === "heatmap" && s.ordering_variable && (
                        // Rotated ordering-variable axis title in the left
                        // margin (R3-Y07, #258). Mockup series-builder.html
                        // :817-822 `.axis-title` — what makes the heatmap read
                        // as a migration map "ordered by X", not stacked rows.
                        <div
                          data-testid="heatmap-axis-title"
                          aria-hidden="true"
                          className="pointer-events-none absolute left-0 top-1/2 z-10 origin-left -translate-y-1/2 -rotate-90 whitespace-nowrap text-[10.5px] uppercase tracking-[0.14em] text-ink-faint"
                        >
                          {s.ordering_variable} &rarr;
                        </div>
                      )}
                      <MultiTracePlot
                        members={members}
                        traces={traces}
                        xDomain={xDomain}
                        onXDomain={setXDomain}
                        groupingMode={groupingMode}
                        sampleIdFor={sampleIdFor}
                        showPeakTicks={showPeakTicks}
                        showPeakLabels={showPeakLabels}
                        xType={scaleMode}
                        workingBandFraction={workingBandFraction}
                        representation={representation}
                        showCrossTraceTracking={trackOn}
                      />
                    </div>
                    <div className="w-[280px] shrink-0" data-testid="series-builder-gutter">
                      <MemberMetaGutter
                        members={members}
                        panelHeight={panelHeight}
                        mode="review"
                        displayLabelByMemberId={displayLabelByMemberId}
                      />
                    </div>
                  </div>

                  {/* Auto figure caption (R8 / B-H). */}
                  <div
                    data-testid="fig-caption"
                    className="mt-3 flex gap-2 border-t border-hair pt-3 text-xs leading-relaxed text-ink-soft"
                  >
                    <span className="font-bold text-ink">Fig.</span>
                    <span>
                      {members.length} 1D integration{members.length === 1 ? "" : "s"}, vertically
                      offset by {offset.toFixed(2)}× the band height
                      {s.ordering_variable ? <>, ordered by {s.ordering_variable}</> : null}. Peak ticks
                      coloured by indexed phase.
                      {/* trackOn caption appendage — mockup `series-builder.html:845` (#251 r1 / N4) */}
                      {trackOn && (
                        <> Tracking lines follow each reflection as it migrates with dose.</>
                      )}
                    </span>
                  </div>
                </>
              )}
            </div>
        </div>
        <SeriesBuilderRail
          collapsed={collapsed}
          onToggleCollapsed={() => setCollapsed((c) => !c)}
          representation={representation}
          onRepresentationChange={setRepresentation}
          orderingVariable={s.ordering_variable}
          offset={offset}
          onOffsetChange={setOffset}
          scaleMode={scaleMode}
          onScaleModeChange={setScaleMode}
          trackOn={trackOn}
          onTrackOnChange={setTrackOn}
          sampleCount={s.samples.length || members.length}
          onAdjustSeries={onStartEdit}
          exportControls={
            <FigureExportControls
              spec={exportSpec}
              filenameStem={exportFilenameStem}
              ariaContext="series figure"
              disabled={exportDisabled}
            />
          }
          {...(editing
            ? { editControls: <SeriesRecipeEditor seriesId={s.id} members={members} /> }
            : {})}
        />
        <OffsetDock
          show={collapsed && representation === "waterfall"}
          value={offset}
          onChange={setOffset}
        />
      </div>
    </ActiveBandProvider>
  );
}
