import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import { PageFrame } from "../components/PageFrame";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSamples,
  useExposures,
  useSetExposureStatus,
  useSelectExposure,
  useAddCorpusSampleTag,
  useRemoveCorpusSampleTag,
} from "../../queries";
import type { Tag } from "../ui";
import { EmptyState, Button } from "../ui";
import { announce } from "../../lib/announce";
import { suppressGlobalKeys } from "../../lib/keys";
import { showToast } from "../../lib/toast";
import { BigFrame } from "../components/BigFrame";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { LoupeSidePanel } from "../components/LoupeSidePanel";
import { PlateHeader } from "../components/PlateHeader";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
  findSampleTagId,
} from "./loupeAdapters";

// Boneyard fixture — a real render with mock props so the headless capture CLI
// measures the greenfield loupe body. image_path:null → DetectorImage takes its
// placeholder branch (a plain rectangle at the frame's fixed aspect).
const FIXTURE_EXPOSURE = {
  id: 0, sample_id: 0, filename: "JC000-001.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null,
  image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
const LOUPE_FIXTURE = (
  <div className="grid grid-cols-[minmax(0,1fr)_286px] gap-7">
    <div className="min-w-0">
      <BigFrame src={null} caption="frame 1 of 1 · kept" />
      <ThumbnailGallery
        exposures={[{ id: 0, src: null, frameNo: 1 }]}
        selectedId={0}
        size="lg"
        align="center"
        className="mt-3"
      />
    </div>
    <LoupeSidePanel
      meta={toMetaEntries(FIXTURE_EXPOSURE, [FIXTURE_EXPOSURE])}
      dropped={false}
      isRepresentative={false}
      tags={[]}
    />
  </div>
);

/**
 * LoupePage (greenfield) — the sample loupe at /samples/loupe/:sampleId.
 * URL-owned: the sample id is the route param, never Zustand `activeSampleId`.
 * Mounts body-only inside the carried CorpusShell <Outlet>.
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  const exposuresQ = useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);
  const isLoading = corpusQ.isLoading || exposuresQ.isLoading;

  // Active exposure — local state, defaulted by defaultExposureId, reset on
  // sample change so the next sample picks its own default.
  const [activeId, setActiveId] = useState<number | undefined>(undefined);
  useEffect(() => { setActiveId(undefined); }, [sampleId]);
  const computedDefault = defaultExposureId(exposures);
  useEffect(() => {
    if (activeId === undefined && computedDefault !== undefined) setActiveId(computedDefault);
  }, [activeId, computedDefault]);
  const activeExposure = exposures.find((e) => e.id === activeId);

  const frameIndex = exposures.findIndex((e) => e.id === activeId);
  const exposurePosition =
    frameIndex >= 0 ? `exposure ${frameIndex + 1} of ${exposures.length}` : "—";

  const setStatus = useSetExposureStatus(hasValidId ? sampleId : 0);
  const setRepresentative = useSelectExposure(hasValidId ? sampleId : 0);
  const addTag = useAddCorpusSampleTag(hasValidId ? sampleId : 0);
  const removeTag = useRemoveCorpusSampleTag(hasValidId ? sampleId : 0);

  const handleDropToggle = useCallback(() => {
    if (!activeExposure) return;
    const dropping = activeExposure.status !== "rejected";
    setStatus.mutate({
      exposureId: activeExposure.id,
      status: dropping ? "rejected" : null,
    });
    // Consequential single-frame status change → visible toast.
    showToast(dropping ? "Frame dropped" : "Frame restored", "success");
  }, [activeExposure, setStatus]);

  const handleSetRepresentative = useCallback(() => {
    if (!activeExposure) return;
    if (activeExposure.selected) {
      // Already the representative — no mutation, no false success toast.
      // Quiet SR-only acknowledgement keeps the R key honest.
      announce("Already the representative frame");
      return;
    }
    setRepresentative.mutate(activeExposure.id);
    showToast("Set as the representative frame", "success");
  }, [activeExposure, setRepresentative]);

  const handleAddTag = useCallback((t: Tag) => {
    addTag.mutate({ key: t.key, value: t.value ?? "" });
  }, [addTag]);

  const handleRemoveTag = useCallback((t: Tag) => {
    if (!sample) return;
    const id = findSampleTagId(sample.tags, t);
    if (id !== undefined) removeTag.mutate(id); // optimistic-add w/o id → no-op (ledger risk)
  }, [removeTag, sample]);

  const flip = useCallback((delta: number) => {
    if (activeId === undefined || exposures.length === 0) return;
    const idx = exposures.findIndex((e) => e.id === activeId);
    if (idx < 0) return;
    const next = Math.min(Math.max(idx + delta, 0), exposures.length - 1);
    setActiveId(exposures[next]!.id);
    // High-frequency, immediately-visible navigation → SR-only (not a toast).
    announce(`Frame ${next + 1} of ${exposures.length}`);
  }, [activeId, exposures]);

  // Thumbnail-click selection routes through the SAME announcement as the
  // keyboard flip so mouse and keyboard navigation are consistent for SR users.
  const selectFrame = useCallback((id: number) => {
    setActiveId(id);
    const idx = exposures.findIndex((e) => e.id === id);
    if (idx >= 0) announce(`Frame ${idx + 1} of ${exposures.length}`);
  }, [exposures]);

  const goBack = useCallback(() => {
    const beamtime = searchParams.get("beamtime");
    navigate(beamtime ? `/samples?beamtime=${beamtime}` : "/samples");
  }, [navigate, searchParams]);

  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      // Modifier chords belong to the browser (⌘R reload, ⌘X cut) — never to
      // the page. The typing/modal suppression is the shared helper.
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      if (suppressGlobalKeys(e)) return;
      if (e.key === "ArrowLeft") flip(-1);
      else if (e.key === "ArrowRight") flip(1);
      else if (e.key === "x" || e.key === "X") handleDropToggle();
      else if (e.key === "r" || e.key === "R") handleSetRepresentative();
      else if (e.key === "Escape") goBack();
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [flip, handleDropToggle, handleSetRepresentative, goBack]);

  if (!corpusQ.isLoading && !sample) {
    return (
      <PageFrame width="loupe" className="px-8 py-7">
        <div data-testid="loupe-page">
          <div data-testid="loupe-not-found">
            <EmptyState
              title="Sample not found"
              body="Nothing in the corpus matches this address."
              action={
                <Button variant="outline" onClick={goBack}>
                  Back to the sheet
                </Button>
              }
            />
          </div>
        </div>
      </PageFrame>
    );
  }

  const isDropped = activeExposure?.status === "rejected";
  // Sample-level truth, not frame-level: the backend's Index-stage resolution
  // never consults status, so a dropped representative still carries forward.
  const representativeDropped = exposures.some(
    (e) => e.selected && e.status === "rejected",
  );

  return (
    <PageFrame width="loupe" className="px-8 py-7">
      <div data-testid="loupe-page">
        <button data-testid="loupe-back" onClick={goBack} className="mb-3.5 text-sm font-semibold text-print-accent hover:underline">
          ← Back to the sheet
        </button>
        <PlateHeader
          title={sample?.display_name ?? sample?.name ?? "—"}
          subtitle={`${sample?.name ?? "—"} · ${exposurePosition}`}
          className="mb-5"
        />
        <Skeleton name="loupe" className="block" loading={isLoading} stagger={50} transition={200}
          fixture={LOUPE_FIXTURE}
          fallback={<div data-testid="loupe-skeleton" className="p-8 text-sm italic text-ink-soft">Loading sample…</div>}>
          <div className="grid grid-cols-[minmax(0,1fr)_286px] gap-7">
            {sample && activeExposure ? (
              <>
                <div className="min-w-0">
                  <BigFrame
                    src={buildExposureImageUrl(activeExposure)}
                    caption={`frame ${frameIndex + 1} of ${exposures.length} · ${isDropped ? "dropped" : "kept"}`}
                    rejected={isDropped}
                  />
                  <ThumbnailGallery
                    exposures={toGalleryExposures(exposures)}
                    {...(activeId !== undefined ? { selectedId: activeId } : {})}
                    onSelect={selectFrame}
                    size="lg"
                    align="center"
                    className="mt-3"
                  />
                </div>
                <LoupeSidePanel
                  meta={toMetaEntries(activeExposure, exposures)}
                  dropped={!!isDropped}
                  isRepresentative={activeExposure.selected}
                  representativeDropped={representativeDropped}
                  tags={toLoupeTags(sample.tags)}
                  onToggleDrop={handleDropToggle}
                  onSetRepresentative={handleSetRepresentative}
                  onAddTag={handleAddTag}
                  onRemoveTag={handleRemoveTag}
                />
              </>
            ) : (
              <div className="col-span-2 p-8 text-sm text-ink-soft">This sample has no exposures.</div>
            )}
          </div>
        </Skeleton>
      </div>
    </PageFrame>
  );
}
