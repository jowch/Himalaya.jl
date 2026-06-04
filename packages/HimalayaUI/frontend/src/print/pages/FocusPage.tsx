import { useMemo, useState } from "react";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { TracePlate } from "../components/TracePlate";
import { DetectorPanel } from "../components/DetectorPanel";
import { CombsPanel, type CombView } from "../components/CombsPanel";
import { AssignmentRail } from "../components/AssignmentRail";
import { AssignmentCart } from "../components/AssignmentCart";
import { PhaseBlock } from "../components/PhaseBlock";
import { CandidateRow, CandidateList } from "../components/CandidateRow";
import { StaleBanner } from "../components/StaleBanner";
import { NotesMargin } from "../components/NotesMargin";
import { CustomIndexModal } from "../components/CustomIndexModal";
import { SpeculativeDialog } from "../components/SpeculativeDialog";
import { ModalShell, Kicker, IconButton, HintText } from "../ui";
import {
  toTraceModel,
  losingPeakIds,
  complementPeakIds,
  toDetectorRings,
  toCombSeries,
  CUSTOM_SYMS,
  customIndexPreview,
} from "./focusAdapters";
import { buildExposureImageUrl } from "./loupeAdapters";
import {
  useCorpusSamples,
  useSamples,
  useExperiment,
  useExposures,
  useTrace,
  usePeaks,
  useIndices,
  useAssignment,
  useUpdateSample,
  useAddPeak,
  useRemovePeak,
  useSetPeakExcluded,
  useAddAssignmentPhase,
  useRemoveAssignmentPhase,
  useCommitCustomIndex,
  useReanalyzeExposure,
  useCreateSpeculative,
  useDeleteIndex,
  useSpeculativeSnap,
} from "../../queries";
import { useAppState } from "../../state";
import { useSyncActiveSampleFromRoute } from "../../hooks/useSyncActiveSampleFromRoute";
import { useAutoPickExposure } from "../../hooks/useAutoPickExposure";
import { useExposureHasPendingPeakOps } from "../../lib/queue/hooks";
import { deriveActiveIndices } from "../../lib/assignment";
import { basisFor } from "../../lib/customIndex";
import { seriesRatio } from "../../lib/seriesRatio";
import { KNOWN_PHASES } from "../../phases";
import type { Trace, IndexEntry } from "../../api";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

/**
 * FocusPage (greenfield) — the Focus workspace at /sample/:sampleId.
 *
 * Mirrors the legacy FocusWorkspaceLayout's 3-column grid (work · rail · notes)
 * but assembled entirely from src/print composites + the carried data layer.
 * The route param seeds Zustand `activeSampleId` (useSyncActiveSampleFromRoute);
 * useAutoPickExposure seeds `activeExposureId`. All hooks are called
 * unconditionally (hooks rule) before any early return.
 */
export function FocusPage(): JSX.Element {
  // ── route → store seeding ──────────────────────────────────────────────────
  useSyncActiveSampleFromRoute();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const notesDrawerOpen = useAppState((s) => s.notesDrawerOpen);
  const closeNotesDrawer = useAppState((s) => s.closeNotesDrawer);
  useAutoPickExposure(activeSampleId);

  // ── data ────────────────────────────────────────────────────────────────────
  const corpusQ = useCorpusSamples();
  const corpusSample =
    activeSampleId !== undefined
      ? corpusQ.data?.find((s) => s.id === activeSampleId)
      : undefined;
  const experimentId = corpusSample?.experiment_id;

  const experimentQ = useExperiment(experimentId ?? 0);
  const samplesQ = useSamples(experimentId ?? 0);
  const exposuresQ = useExposures(activeSampleId);

  const traceQ = useTrace(activeExposureId);
  const peaksQ = usePeaks(activeExposureId);
  const indicesQ = useIndices(activeExposureId);
  const assignmentQ = useAssignment(activeExposureId);

  // CACHE-COHERENCE: the notes textarea must read from the SAME cache the
  // update mutator patches — the experiment-scoped samples list, NOT the corpus
  // list (mirrors the legacy FocusWorkspaceLayout note).
  const notesSample =
    activeSampleId !== undefined && experimentId !== undefined
      ? samplesQ.data?.find((s) => s.id === activeSampleId)
      : undefined;

  // ── mutators (all exposure-scoped, except the sample-scoped notes save) ──────
  const addPeak = useAddPeak(activeExposureId ?? 0);
  const removePeak = useRemovePeak(activeExposureId ?? 0);
  const setPeakExcluded = useSetPeakExcluded(activeExposureId ?? 0);
  const addAssignmentPhase = useAddAssignmentPhase(activeExposureId ?? 0);
  const removeAssignmentPhase = useRemoveAssignmentPhase(activeExposureId ?? 0);
  const commitCustomIndex = useCommitCustomIndex(activeExposureId ?? 0);
  const reanalyze = useReanalyzeExposure(activeExposureId ?? 0);
  const createSpeculative = useCreateSpeculative(activeExposureId ?? 0);
  const deleteIndex = useDeleteIndex(activeExposureId ?? 0);
  const updateSample = useUpdateSample(experimentId ?? 0, activeSampleId ?? 0);
  const pendingPeakOps = useExposureHasPendingPeakOps(activeExposureId);

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

  // speculative dialog
  const [specOpen, setSpecOpen] = useState(false);
  const [specPhase, setSpecPhase] = useState<string>("Lamellar");
  const [anchorPeakId, setAnchorPeakId] = useState<number | undefined>(undefined);
  const [anchorRatio, setAnchorRatio] = useState<number>(1);
  const [included, setIncluded] = useState<Record<number, boolean>>({});

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

  // Stale-index count: indices whose snapshot input-hash differs from the
  // active exposure's analysis hash (the StaleIndicesBanner contract).
  const expectedHash = activeExposure?.analysis_inputs_hash ?? null;
  const staleCount = expectedHash
    ? indices.filter((i) => i.inputs_hash !== expectedHash).length
    : 0;

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

  // speculative snap (gated by the carried hook)
  const snapQ = useSpeculativeSnap(
    activeExposureId,
    specPhase,
    anchorPeakId,
    anchorRatio,
  );
  const snap = snapQ.data ?? [];

  const isLoading =
    corpusQ.isLoading ||
    (activeExposureId !== undefined && (traceQ.isLoading || peaksQ.isLoading));

  // ── early states ─────────────────────────────────────────────────────────────
  if (!corpusQ.isLoading && !corpusSample) {
    return (
      <PageFrame width="focus" className="px-8 py-7">
        <div
          data-testid="focus-not-found"
          className="rounded border border-hair-strong p-8 text-sm text-ink-faint"
        >
          Sample not found.
        </div>
      </PageFrame>
    );
  }

  const noExposure =
    !exposuresQ.isLoading && activeExposureId === undefined;

  // ── speculative helpers ──────────────────────────────────────────────────────
  function commitCustom(): void {
    commitCustomIndex.mutate(customSym, basisFor(customSym, Number(customParam)));
    setCustomOpen(false);
  }

  function createSpec(): void {
    if (anchorPeakId === undefined) return;
    const additional = snap
      .filter(
        (s) =>
          !s.is_anchor &&
          s.suggested_peak_id !== null &&
          included[s.ratio_position],
      )
      .map((s) => ({
        ratio_position: s.ratio_position,
        peak_id: s.suggested_peak_id!,
      }));
    createSpeculative.mutate({
      phase: specPhase,
      anchor_peak_id: anchorPeakId,
      anchor_ratio: anchorRatio,
      additional,
      active: false,
    });
    setSpecOpen(false);
  }

  // ── rail building ────────────────────────────────────────────────────────────
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
          onRemove={() => removeAssignmentPhase.mutate(ix.id)}
        />
      ))}
    </AssignmentCart>
  );

  function candidateRow(ix: IndexEntry, deletable = false): JSX.Element {
    const selected = memberIds.has(ix.id);
    return (
      // Placement-only wrapper so the candidate-hover preview (losing-peak dim)
      // can hang off mouseEnter/leave — CandidateRow itself is a presentational
      // primitive with no hover hook.
      <div
        key={ix.id}
        className="flex items-center gap-1.5"
        onMouseEnter={() => setPreviewIndexId(ix.id)}
        onMouseLeave={() => setPreviewIndexId(undefined)}
      >
        <CandidateRow
          className="flex-1 min-w-0"
          phase={ix.phase}
          score={ix.score}
          why={`explains ${ix.peaks.length} peaks${selected ? " · in the call" : ""}`}
          selected={selected}
          {...(ix.bonnet?.consistent ? { bonnet: true } : {})}
          onToggle={() =>
            selected
              ? removeAssignmentPhase.mutate(ix.id)
              : addAssignmentPhase.mutate(ix.id)
          }
        />
        {deletable && (
          <IconButton
            label={`Delete speculative index ${ix.id}`}
            tone="danger"
            data-testid={`spec-delete-${ix.id}`}
            onClick={() => deleteIndex.mutate(ix.id)}
            className="shrink-0"
          >
            <svg viewBox="0 0 16 16" className="h-3.5 w-3.5" aria-hidden="true">
              <path
                d="M3.5 4.5h9M6 4.5V3.2a1 1 0 011-1h2a1 1 0 011 1v1.3M5 4.5l.55 8a1 1 0 001 .93h2.9a1 1 0 001-.93l.55-8"
                fill="none"
                stroke="currentColor"
                strokeWidth="1.3"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </IconButton>
        )}
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
      {speculatives.map((ix) => candidateRow(ix, true))}
      <button
        type="button"
        data-testid="add-speculative-button"
        onClick={() => setSpecOpen(true)}
        className="mt-1 w-full rounded-md border border-dashed border-hair py-1.5 text-caption text-ink-faint transition-colors hover:bg-paper-sunk hover:text-ink"
      >
        + Add speculative
      </button>
    </CandidateList>
  );

  // ── notes (shared by the xl margin + the < xl drawer) ────────────────────────
  const notesBody =
    notesSample !== undefined ? (
      <NotesMargin
        notes={notesSample.notes ?? null}
        onSaveNotes={(notes) => updateSample.mutate({ notes })}
        onHoverQ={setHoveredQ}
      />
    ) : null;

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
      <SpeculativeDialog
        open={specOpen}
        onClose={() => setSpecOpen(false)}
        phases={KNOWN_PHASES as readonly string[]}
        phase={specPhase}
        onPhaseChange={setSpecPhase}
        peaks={peaks.map((p) => ({ id: p.id, q: p.q, source: p.source }))}
        anchorPeakId={anchorPeakId}
        onAnchorChange={setAnchorPeakId}
        anchorRatio={anchorRatio}
        onAnchorRatioChange={setAnchorRatio}
        snap={snap}
        included={included}
        onToggleIncluded={(rp) =>
          setIncluded((prev) => ({ ...prev, [rp]: !prev[rp] }))
        }
        snapLoading={snapQ.isLoading}
        blocked={pendingPeakOps}
        saving={createSpeculative.isPending}
        onCreate={createSpec}
      />
    </>
  );

  // ── layout ───────────────────────────────────────────────────────────────────
  return (
    <PageFrame width="focus" className="px-8 py-7">
      <Skeleton
        name="focus"
        className="block"
        loading={isLoading}
        stagger={50}
        transition={200}
        fallback={
          <div
            data-testid="focus-skeleton"
            className="p-8 text-sm italic text-ink-faint"
          >
            Loading workspace…
          </div>
        }
      >
        <div
          data-testid="focus-workspace"
          className="grid grid-cols-1 gap-5 xl:grid-cols-[1fr_348px_250px]"
        >
          {/* work area */}
          <div className="flex min-h-0 flex-col gap-5">
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
                onAddPeak: (q) => addPeak.mutate(q),
                onClickPeak: (id, alt) =>
                  alt
                    ? setPeakExcluded.mutate({ peakId: id, excluded: true })
                    : removePeak.mutate(id),
                onReset: () => setXDomain(null),
              }}
            />
            <div className="grid grid-cols-1 gap-5 lg:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
              <DetectorPanel
                src={detectorSrc}
                rings={rings}
                {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                onHoverQ={setHoveredQ}
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

          {/* rail */}
          <div className="flex min-h-0 flex-col gap-4">
            {staleCount > 0 && (
              <StaleBanner
                staleCount={staleCount}
                pending={pendingPeakOps}
                onReanalyze={() => reanalyze.mutate({})}
              />
            )}
            {noExposure ? (
              <div className="p-8 text-sm text-ink-faint">
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

          {/* notes margin (xl+) */}
          {notesBody && <div className="hidden xl:block">{notesBody}</div>}
        </div>
      </Skeleton>

      {/* notes drawer (< xl) */}
      {notesBody && (
        <div className="xl:hidden">
          <ModalShell
            open={notesDrawerOpen}
            onClose={closeNotesDrawer}
            variant="drawer"
            testId="focus-notes-drawer"
            aria-label="Notes"
          >
            <div className="flex items-center justify-between border-b border-hair px-4 py-2">
              <Kicker as="span" tone="faint">
                Notes
              </Kicker>
              <IconButton
                label="Close notes"
                dismiss
                tone="ghost"
                onClick={closeNotesDrawer}
              />
            </div>
            {notesBody}
          </ModalShell>
        </div>
      )}

      {modals}
    </PageFrame>
  );
}
