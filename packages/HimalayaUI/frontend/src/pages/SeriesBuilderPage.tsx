import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useSeries, useMemberTraces, useMemberTracesLoading,
  useMemberExposures, useMemberSamples,
} from "../queries";
import { MultiTracePlot } from "../components/MultiTracePlot";
import { MemberMetaGutter } from "../components/MemberMetaGutter";
import { GroupingModeToggle } from "../components/GroupingModeToggle";
import { AnnotationToggles } from "../components/AnnotationToggles";
import { ActiveBandProvider } from "../components/ActiveBandContext";
import { FigureExportControls } from "../components/FigureExportControls";
import { SeriesBuilderRail } from "../components/SeriesBuilderRail";
import { SeriesRecipeEditor } from "../components/SeriesRecipeEditor";
import { HintText } from "../components/ui";
import type { Representation } from "../components/RepresentationToggle";
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
  const title = s && s.title.trim() !== "" ? s.title : "Untitled series";

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
        {s && (
          <>
            <header
              data-testid="series-builder-header"
              data-editing={String(editing)}
              className="shrink-0 flex items-start justify-between px-6 pt-5"
            >
              <div>
                <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
                  Series
                </div>
                <h1 className="font-medium text-ink">{title}</h1>
              </div>
              {editing ? (
                <span data-testid="series-builder-editing-badge" className="text-xs text-ink-faint">
                  Editing
                </span>
              ) : (
                <button
                  type="button"
                  data-testid="series-builder-edit"
                  onClick={() => startSeriesDraft(s)}
                  className="rounded border border-hair px-2 py-1 text-sm text-ink hover:bg-paper-sunk"
                >
                  Edit
                </button>
              )}
            </header>
            {s.members.length === 0 ? (
              <div
                data-testid="series-builder-empty"
                className="flex-1 grid place-items-center text-sm text-ink-faint"
              >
                This series has no members yet.
              </div>
            ) : (
              <SeriesBuilderBody series={s} editing={editing} />
            )}
          </>
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
  { series: s, editing }: { series: Series; editing: boolean },
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
  // onXDomain. Compare keeps this in Zustand `compareXDomains[comparisonId]`,
  // but that Record is keyed by COMPARISON id — reusing it for a series id
  // would namespace-collide. A read surface's pan/zoom is local, so hold it
  // in local useState. null = full data range.
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
  // (waterfall live; heatmap deferred to #208).
  const [collapsed, setCollapsed] = useState(false);
  const [representation, setRepresentation] = useState<Representation>("waterfall");

  // Annotation toggles live in Zustand (shared with AnnotationToggles); the
  // export spec must reflect their current value so the figure matches the
  // on-screen plot.
  const showPeakTicks = useAppState((st) => st.showPeakTicks);
  const showPeakLabels = useAppState((st) => st.showPeakLabels);

  // Figure export — mirror Compare's spec thunk (evaluated at click time so
  // it captures fresh xDomain / toggles). No experiment scope on this surface,
  // so experimentName is omitted and the filename stem uses the series title.
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
  }), [
    members, traces, s.title, xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
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
        <div className="flex-1 min-h-0 flex flex-col p-4 gap-3" data-testid="series-builder-plot">
          <div className="flex items-center gap-3" data-testid="series-builder-controls">
            <GroupingModeToggle mode={groupingMode} onChange={setGroupingMode} />
            <AnnotationToggles />
          </div>
          <div className="flex-1 min-h-0 flex flex-row gap-2">
            <div ref={plotColRef} className="flex-1 min-w-0">
              <MultiTracePlot
                members={members}
                traces={traces}
                xDomain={xDomain}
                onXDomain={setXDomain}
                groupingMode={groupingMode}
                sampleIdFor={sampleIdFor}
                showPeakTicks={showPeakTicks}
                showPeakLabels={showPeakLabels}
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
        </div>
        <SeriesBuilderRail
          collapsed={collapsed}
          onToggleCollapsed={() => setCollapsed((c) => !c)}
          representation={representation}
          onRepresentationChange={setRepresentation}
          orderingVariable={s.ordering_variable}
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
      </div>
    </ActiveBandProvider>
  );
}
