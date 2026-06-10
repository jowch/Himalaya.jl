import { useMemo, useState, useCallback } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { TracePlate } from "../components/TracePlate";
import { DetectorPanel } from "../components/DetectorPanel";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { CombsPanel, type CombView } from "../components/CombsPanel";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { AssignmentRail } from "../components/AssignmentRail";
import { AssignmentCart } from "../components/AssignmentCart";
import { PhaseBlock } from "../components/PhaseBlock";
import { CandidateRow, CandidateList } from "../components/CandidateRow";
import { CustomIndexModal } from "../components/CustomIndexModal";
import { HintText, EmptyState, Button } from "../ui";
import { ExportButton } from "../components/ExportButton";
import { useFigureExport } from "../components/useFigureExport";
import { buildTraceExportSpec } from "../../lib/figure-export/adapters/traceAdapter";
import type { ExportSpec } from "../../lib/figure-export/types";
import {
  toTraceModel,
  losingPeakIds,
  complementPeakIds,
  toDetectorRings,
  toCombSeries,
  CUSTOM_SYMS,
  customIndexPreview,
} from "./focusAdapters";
import { buildExposureImageUrl, toGalleryExposures } from "./loupeAdapters";
import {
  useCorpusSamples,
  useExperiment,
  useExposures,
  useTrace,
  usePeaks,
  useIndices,
  useAssignment,
  useAddPeak,
  useRemovePeak,
  useSetPeakExcluded,
  useAddAssignmentPhase,
  useRemoveAssignmentPhase,
  useCommitCustomIndex,
} from "../../queries";
import { useAppState } from "../../state";
import { useSyncActiveSampleFromRoute } from "../../hooks/useSyncActiveSampleFromRoute";
import { useAutoPickExposure } from "../../hooks/useAutoPickExposure";
import { deriveActiveIndices } from "../../lib/assignment";
import { basisFor } from "../../lib/customIndex";
import { seriesRatio } from "../../lib/seriesRatio";
import { announce } from "../../lib/announce";
import { showToast } from "../../lib/toast";
import type { Trace, IndexEntry } from "../../api";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

// Boneyard fixture — a real render with mock props so the headless capture CLI
// measures the greenfield Focus body. Uses static composites with no-op handlers
// and a small synthetic trace so the plot draws real geometry.
const FIXTURE_TRACE = {
  trace: {
    q: [0.02, 0.04, 0.06, 0.1, 0.15, 0.2],
    I: [12000, 8000, 4800, 2400, 1200, 600],
    sigma: [80, 55, 35, 20, 12, 7],
  },
  peaks: [] as PlotPeak[],
  phase: null,
};
const FOCUS_FIXTURE = (
  <div className="grid grid-cols-1 lg:grid-cols-[minmax(0,1fr)_350px]">
    {/* work column */}
    <div className="min-w-0 px-8 pt-7 pb-13">
      <div className="mx-auto flex min-w-0 max-w-[1180px] flex-col gap-[18px]">
        <TracePlate
          kicker="Integration"
          title="Sample"
          subtitle="JC000 · Experiment · frame-001"
          trace={FIXTURE_TRACE}
          scale="log"
          onScaleChange={() => {}}
          onAutoFit={() => {}}
          onToggleAddPeak={() => {}}
        />
        <div className="grid grid-cols-1 gap-5 lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
          <DetectorPanel src={null} />
          <div className="hidden lg:flex min-h-0 flex-col">
            <CombsPanel
              assigned={[]}
              view="comb"
              onViewChange={() => {}}
            />
          </div>
        </div>
      </div>
    </div>

    {/* rail column = AssignmentRail (owns the rail shell) */}
    <AssignmentRail
      assignment={
        <AssignmentCart>
          <PhaseBlock phase="Pn3m" score={0.87} meta="a = 142 Å · 5 reflections" />
        </AssignmentCart>
      }
      candidates={
        <CandidateList>
          <CandidateRow phase="Pn3m" score={0.87} why="explains 5 peaks · in the call" selected />
          <CandidateRow phase="Im3m" score={0.61} why="explains 3 peaks" />
        </CandidateList>
      }
      candidatesNote="Check every phase that is present; a sample can hold more than one."
    />
  </div>
);

/**
 * FocusPage (greenfield) — the Focus workspace at /sample/:sampleId.
 *
 * Mirrors the 2-column workspace grid (work · rail)
 * but assembled entirely from src/print composites + the carried data layer.
 * The route param seeds Zustand `activeSampleId` (useSyncActiveSampleFromRoute);
 * useAutoPickExposure seeds `activeExposureId`. All hooks are called
 * unconditionally (hooks rule) before any early return.
 */
export function FocusPage(): JSX.Element {
  const navigate = useNavigate();

  // ── route → store seeding ──────────────────────────────────────────────────
  useSyncActiveSampleFromRoute();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  useAutoPickExposure(activeSampleId);

  // ── data ────────────────────────────────────────────────────────────────────
  const corpusQ = useCorpusSamples();
  const corpusSample =
    activeSampleId !== undefined
      ? corpusQ.data?.find((s) => s.id === activeSampleId)
      : undefined;
  const experimentId = corpusSample?.experiment_id;

  const experimentQ = useExperiment(experimentId ?? 0);
  const exposuresQ = useExposures(activeSampleId);

  const traceQ = useTrace(activeExposureId);
  const peaksQ = usePeaks(activeExposureId);
  const indicesQ = useIndices(activeExposureId);
  const assignmentQ = useAssignment(activeExposureId);

  // ── mutators (all exposure-scoped) ───────────────────────────────────────────
  const addPeak = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);
  const setPeakExcluded = useSetPeakExcluded(activeExposureId ?? 0);
  const addAssignmentPhase = useAddAssignmentPhase(activeExposureId ?? 0);
  const removeAssignmentPhase = useRemoveAssignmentPhase(activeExposureId ?? 0);
  const commitCustomIndex = useCommitCustomIndex(activeExposureId ?? 0);

  // ── page-owned state ─────────────────────────────────────────────────────────
  const [scale, setScale] = useState<"log" | "lin">("log");
  const [addArmed, setAddArmed] = useState(false);
  const [xDomain, setXDomain] = useState<[number, number] | null>(null);
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  const [combView, setCombView] = useState<CombView>("comb");
  const [previewIndexId, setPreviewIndexId] = useState<number | undefined>(undefined);

  // custom-index modal
  const [customOpen, setCustomOpen] = useState(false);
  const [customSym, setCustomSym] = useState<string>(CUSTOM_SYMS[0]!.name);
  const [customParam, setCustomParam] = useState<string>(
    String(CUSTOM_SYMS[0]!.min),
  );

  // ── derived ──────────────────────────────────────────────────────────────────
  const peaks = peaksQ.data ?? [];
  const indices = indicesQ.data ?? [];
  const assignment = assignmentQ.data;
  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);
  const activeExposure = exposures.find((e) => e.id === activeExposureId);

  const activeIndices = useMemo(
    () => deriveActiveIndices(assignment, indices),
    [assignment, indices],
  );
  const memberIds = useMemo(
    () => new Set(activeIndices.map((ix) => ix.id)),
    [activeIndices],
  );

  // Primary assigned phase: the first active index when the assignment is
  // `indexed` (mirrors the legacy trace-colour derivation, which reads the
  // active set off the assignment cart). Null otherwise → ink-faint trace.
  const phase =
    assignment?.state === "indexed" && activeIndices.length > 0
      ? activeIndices[0]!.phase
      : null;

  const traceModel = useMemo(
    () => toTraceModel(traceQ.data ?? EMPTY_TRACE, peaks, phase),
    [traceQ.data, peaks, phase],
  );

  // Candidate-hover preview → losing-peak dim. The hovered candidate's claimed
  // peaks would orphan some active-phase peaks; those "losing" peaks dim while
  // the rest stay highlighted (TracePlot KEEPS the complement set).
  const previewIx = indices.find((i) => i.id === previewIndexId);
  const highlight = useMemo(() => {
    if (!previewIx) return undefined;
    const losing = losingPeakIds(previewIx, activeIndices);
    if (losing.size === 0) return undefined;
    return complementPeakIds(
      peaks.map((p) => p.id),
      losing,
    );
  }, [previewIx, activeIndices, peaks]);

  const rings = useMemo(
    () => toDetectorRings(activeIndices, peaks),
    [activeIndices, peaks],
  );
  const { assigned: combAssigned, leftover: combLeftover } = useMemo(
    () => toCombSeries(activeIndices, peaks),
    [activeIndices, peaks],
  );

  // Detector image url for the active exposure (loupe's builder → cache-coherent).
  const detectorSrc = activeExposure
    ? buildExposureImageUrl(activeExposure)
    : null;

  const exposureLabel = activeExposure?.filename
    ? activeExposure.filename.replace(/\.[^.]+$/, "")
    : null;
  const sampleName =
    corpusSample?.display_name ?? corpusSample?.name ?? "—";
  const subtitle = [
    corpusSample?.name,
    experimentQ.data?.name,
    exposureLabel,
  ]
    .filter(Boolean)
    .join(" · ");

  // custom-index live preview
  const customMeta = CUSTOM_SYMS.find((s) => s.name === customSym) ?? CUSTOM_SYMS[0]!;
  const observedQs = peaks.filter((p) => !p.excluded).map((p) => p.q);
  const { previewSeries, fit } = customIndexPreview(
    customSym,
    Number(customParam),
    observedQs,
  );

  const isLoading =
    corpusQ.isLoading ||
    (activeExposureId !== undefined && (traceQ.isLoading || peaksQ.isLoading));

  // ── figure export ─────────────────────────────────────────────────────────────
  const exportSpec = useCallback((): ExportSpec => buildTraceExportSpec({
    trace: traceQ.data ?? EMPTY_TRACE,
    peaks,
    activeGroupIndices: activeIndices,
    experimentName: experimentQ.data?.name ?? "",
    sampleName,
    exposureLabel: exposureLabel ?? "",
    xDomain,
    yDomain: null,
    xType: scale === "log" ? "log" : "linear",
    ...(experimentQ.data?.q_units ? { qUnits: experimentQ.data.q_units } : {}),
  }), [traceQ.data, peaks, activeIndices, experimentQ.data, sampleName, exposureLabel, xDomain, scale]);

  const filenameStem = `${sampleName} ${exposureLabel ?? ""}`.trim();
  const fx = useFigureExport(exportSpec, filenameStem, "trace plot");

  // ── early states ─────────────────────────────────────────────────────────────
  if (!corpusQ.isLoading && !corpusSample) {
    return (
      <div className="px-8 py-7">
        <div data-testid="focus-not-found">
          <EmptyState
            title="Sample not found"
            body="Nothing in the corpus matches this address."
            action={
              <Button variant="outline" onClick={() => navigate("/samples")}>
                Back to the contact sheet
              </Button>
            }
          />
        </div>
      </div>
    );
  }

  const noExposure =
    !exposuresQ.isLoading && activeExposureId === undefined;

  // ── custom-index helper ──────────────────────────────────────────────────────
  function commitCustom(): void {
    commitCustomIndex.mutate(customSym, basisFor(customSym, Number(customParam)));
    setCustomOpen(false);
    // Consequential: a custom hypothesis was added to the call → visible toast.
    showToast(`${customSym} index added`, "success");
  }

  // ── rail building ────────────────────────────────────────────────────────────
  // Speculative-kind indices are the output of the custom-index builder (it
  // persists a speculative index + adds it to the assignment); they render as
  // ordinary toggle candidate rows alongside the auto candidates — the mockup's
  // candidate list is uniform/toggle-only (no per-row delete). The custom-index
  // modal is the hypothesis tool; the cart's remove (PhaseBlock onRemove) takes
  // a phase back out of the call.
  const speculatives = indices.filter((i) => i.kind === "speculative");
  const candidatePool = indices.filter((i) => i.kind !== "speculative");

  function phaseMeta(ix: IndexEntry): string {
    const lattice =
      ix.lattice_d != null ? `a = ${ix.lattice_d.toFixed(0)} Å · ` : "";
    return `${lattice}${ix.peaks.length} reflections`;
  }

  const cart = (
    <AssignmentCart onCustomIndex={() => setCustomOpen(true)}>
      {activeIndices.map((ix) => (
        <PhaseBlock
          key={ix.id}
          phase={ix.phase}
          score={ix.score ?? 0}
          meta={phaseMeta(ix)}
          series={seriesRatio(ix.phase, ix.peaks.map((p) => p.ratio_position)) || undefined}
          onRemove={() => {
            removeAssignmentPhase.mutate(ix.id);
            announce(`${ix.phase} removed from the call`);
          }}
        />
      ))}
    </AssignmentCart>
  );

  function candidateRow(ix: IndexEntry): JSX.Element {
    const selected = memberIds.has(ix.id);
    return (
      // Placement-only wrapper so the candidate-hover preview (losing-peak dim)
      // can hang off mouseEnter/leave — CandidateRow itself is a presentational
      // primitive with no hover hook.
      <div
        key={ix.id}
        onMouseEnter={() => setPreviewIndexId(ix.id)}
        onMouseLeave={() => setPreviewIndexId(undefined)}
      >
        <CandidateRow
          phase={ix.phase}
          score={ix.score}
          why={`explains ${ix.peaks.length} peaks${selected ? " · in the call" : ""}`}
          selected={selected}
          {...(ix.bonnet?.consistent ? { bonnet: true } : {})}
          onToggle={() => {
            if (selected) {
              removeAssignmentPhase.mutate(ix.id);
              announce(`${ix.phase} removed from the call`);
            } else {
              addAssignmentPhase.mutate(ix.id);
              announce(`${ix.phase} added to the call`);
            }
          }}
        />
      </div>
    );
  }

  const candidates = (
    <CandidateList>
      {candidatePool.length === 0 ? (
        <HintText>No candidate indexings.</HintText>
      ) : (
        candidatePool.map((ix) => candidateRow(ix))
      )}
      {speculatives.map((ix) => candidateRow(ix))}
    </CandidateList>
  );

  // ── modals (mounted regardless of layout branch) ─────────────────────────────
  const modals = (
    <>
      <CustomIndexModal
        open={customOpen}
        symmetries={CUSTOM_SYMS.map((s) => s.name)}
        symmetry={customSym}
        onSymmetryChange={setCustomSym}
        paramName={customMeta.paramName}
        paramMin={customMeta.min}
        paramMax={customMeta.max}
        {...(customMeta.step !== undefined ? { paramStep: customMeta.step } : {})}
        unit={customMeta.unit}
        paramValue={customParam}
        onParamChange={setCustomParam}
        previewSeries={previewSeries}
        observed={observedQs}
        fit={fit}
        onAdd={commitCustom}
        onCancel={() => setCustomOpen(false)}
        onClose={() => setCustomOpen(false)}
      />
    </>
  );

  // ── layout ───────────────────────────────────────────────────────────────────
  // Full-bleed per the focus-plot mockup: a [work 1fr · rail 350px] grid with
  // NO outer max-width — the work column centres its own content at 1180px and
  // the rail (bg-paper-sunk + left hairline) pins to the right edge.
  return (
    <>
      <Skeleton
        name="focus"
        className="block"
        loading={isLoading}
        stagger={50}
        transition={200}
        fixture={FOCUS_FIXTURE}
        fallback={
          <div
            data-testid="focus-skeleton"
            className="p-8 text-sm italic text-ink-soft"
          >
            Loading workspace…
          </div>
        }
      >
        <div
          data-testid="focus-workspace"
          className="grid grid-cols-1 lg:grid-cols-[minmax(0,1fr)_350px]"
        >
          {/* work column — full-bleed; inner content capped at 1180px (mockup .work / .work-inner) */}
          <div className="min-w-0 px-8 pt-7 pb-13">
           <div className="mx-auto flex min-w-0 max-w-[1180px] flex-col gap-[18px]">
            <TracePlate
              kicker="Integration"
              title={sampleName}
              subtitle={subtitle}
              trace={traceModel}
              scale={scale}
              onScaleChange={setScale}
              addPeakArmed={addArmed}
              onToggleAddPeak={() => setAddArmed((v) => !v)}
              onAutoFit={() => setXDomain(null)}
              xDomain={xDomain}
              {...(hoveredQ !== undefined ? { hoveredQ } : {})}
              onHoverQ={setHoveredQ}
              {...(highlight !== undefined ? { highlightPeakIds: highlight } : {})}
              interaction={{
                onXDomain: setXDomain,
                onAddPeak: (q) => {
                  addPeak.mutate(q);
                  announce("Peak added");
                },
                onClickPeak: (id, alt) => {
                  if (alt) {
                    setPeakExcluded.mutate({ peakId: id, excluded: true });
                    announce("Peak excluded");
                  } else {
                    removePeak.mutate(id);
                    announce("Peak removed");
                  }
                },
                onReset: () => setXDomain(null),
              }}
              actions={
                <ExportButton
                  {...fx}
                  ariaContext="trace plot"
                  disabled={!traceQ.data || peaks.length === 0}
                />
              }
            />
            <div className="grid grid-cols-1 gap-5 lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
              <DetectorPanel
                src={detectorSrc}
                rings={rings}
                {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                onHoverQ={setHoveredQ}
                tools={
                  exposures.length > 0 ? (
                    <ThumbnailGallery
                      exposures={toGalleryExposures(exposures)}
                      {...(activeExposureId !== undefined ? { selectedId: activeExposureId } : {})}
                      onSelect={setActiveExposure}
                      size="xs"
                    />
                  ) : undefined
                }
              />
              <div className="hidden lg:flex min-h-0 flex-col">
                <CombsPanel
                  assigned={combAssigned}
                  leftover={combLeftover}
                  view={combView}
                  onViewChange={setCombView}
                  {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                  onHoverQ={setHoveredQ}
                />
              </div>
            </div>
           </div>
          </div>

          {/* rail column = AssignmentRail itself (it owns the rail shell:
              bg-paper-sunk + left hairline + padding + scroll). The empty state
              mirrors that shell so the column looks identical with no exposures. */}
          {noExposure ? (
            <div className="flex min-h-0 flex-col border-l border-hair bg-paper-sunk p-5 text-sm text-ink-soft">
              This sample has no exposures.
            </div>
          ) : (
            <AssignmentRail
              assignmentCount={activeIndices.length || undefined}
              assignment={cart}
              candidates={candidates}
              candidatesNote="Check every phase that is present; a sample can hold more than one. Candidates that explain the same peaks swap; independent phases coexist."
            />
          )}
        </div>
      </Skeleton>

      {modals}
    </>
  );
}
