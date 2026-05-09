import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import {
  useExperiment,
  useExposures,
  useSamples,
  useSetExposureStatus,
  useSelectExposure,
  useAddExposureTag,
  useRemoveExposureTag,
  useAddSampleTag,
  useRemoveSampleTag,
  useUpdateSample,
} from "../queries";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { DetectorImageCard } from "../components/DetectorImageCard";
import { SampleMetadataCard } from "../components/SampleMetadataCard";
import { ChatCard } from "../components/ChatCard";
import { WorkspaceGrid } from "../components/WorkspaceGrid";
import { WarmAddMenu } from "../components/WarmAddMenu";
import type { Exposure, Sample } from "../api";

const DETECTOR_IMAGE_FIXTURE_EXPOSURE: Exposure = {
  id: 0,
  sample_id: 0,
  filename: "JC001-001.dat",
  kind: "file",
  selected: false,
  status: "accepted",
  image_path: null,
  image_version: "",
  tags: [],
  sources: [],
  trace_hash: null,
  analysis_inputs_hash: null,
};

const DETECTOR_IMAGE_FIXTURE = (
  <DetectorImageCard
    exposure={DETECTOR_IMAGE_FIXTURE_EXPOSURE}
    onSetStatus={() => {}}
    onSetIndexing={() => {}}
    onAddTag={() => {}}
    onRemoveTag={() => {}}
  />
);

const THUMBNAIL_GALLERY_FIXTURE_EXPOSURES: Exposure[] = Array.from(
  { length: 6 },
  (_, i) => ({
    id: i + 1,
    sample_id: 0,
    filename: `JC001-00${i + 1}.dat`,
    kind: "file" as const,
    selected: i === 0,
    status: "accepted" as const,
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
  }),
);

const THUMBNAIL_GALLERY_FIXTURE = (
  <ThumbnailGallery
    exposures={THUMBNAIL_GALLERY_FIXTURE_EXPOSURES}
    selectedId={1}
    onSelect={() => {}}
    className="h-full"
  />
);

const SAMPLE_METADATA_FIXTURE_SAMPLE: Sample = {
  id: 0,
  experiment_id: 0,
  name: "JC001",
  display_name: "DOPE 70%",
  notes: null,
  tags: [{ id: 1, key: "lipid", value: "DOPE", source: "manifest" }],
};

const SAMPLE_METADATA_FIXTURE = (
  <SampleMetadataCard
    sample={SAMPLE_METADATA_FIXTURE_SAMPLE}
    experimentName="Experiment A"
    exposureSummary={{ total: 8, accepted: 6, rejected: 1 }}
    onUpdateSample={() => {}}
    onAddTag={() => {}}
    onRemoveTag={() => {}}
  />
);

export function InspectPage(): JSX.Element {
  const username     = useAppState((s) => s.username);
  const experimentId = useAppState((s) => s.activeExperimentId);
  const sampleId     = useAppState((s) => s.activeSampleId);
  const openModal    = useAppState((s) => s.openNavModal);

  // Auto-open nav modal if no sample selected (mirrors IndexPage behaviour)
  const autoOpenedRef = useRef(false);
  useEffect(() => {
    if (autoOpenedRef.current) return;
    if (username === undefined) return;
    if (experimentId === undefined) {
      autoOpenedRef.current = true;
      openModal("experiment");
    } else if (sampleId === undefined) {
      autoOpenedRef.current = true;
      openModal("sample");
    }
  }, [username, experimentId, sampleId, openModal]);

  const experimentQ = useExperiment(experimentId ?? 0);
  const exposuresQ  = useExposures(sampleId);
  const samplesQ    = useSamples(experimentId ?? 0);
  const exposures   = exposuresQ.data ?? [];
  const sample      = samplesQ.data?.find((s) => s.id === sampleId);
  const experimentName =
    experimentQ.data?.name ?? experimentQ.data?.path ?? undefined;

  // Default: indexing-marked → first accepted → first
  const defaultId = useMemo(() => {
    const indexing = exposures.find((e) => e.selected);
    if (indexing) return indexing.id;
    const firstAccepted = exposures.find((e) => e.status === "accepted");
    if (firstAccepted) return firstAccepted.id;
    return exposures[0]?.id;
  }, [exposures]);

  const [viewingId, setViewingId] = useState<number | undefined>(undefined);

  // Reset on sample switch so the next render picks the new sample's default
  // (indexing → first accepted → first). Without this, viewingId stays pinned
  // to an exposure id from the previous sample.
  useEffect(() => {
    setViewingId(undefined);
  }, [sampleId]);

  useEffect(() => {
    if (viewingId === undefined && defaultId !== undefined) {
      setViewingId(defaultId);
    }
  }, [defaultId, viewingId]);

  const viewingExposure = exposures.find((e) => e.id === viewingId);

  const setStatus    = useSetExposureStatus(sampleId ?? 0);
  const setIndexing  = useSelectExposure(sampleId ?? 0);
  const addExpTag    = useAddExposureTag(sampleId ?? 0, viewingId ?? 0);
  const rmExpTag     = useRemoveExposureTag(sampleId ?? 0, viewingId ?? 0);
  const updateSample = useUpdateSample(experimentId ?? 0, sampleId ?? 0);
  const addSampleTag = useAddSampleTag(experimentId ?? 0, sampleId ?? 0);
  const rmSampleTag  = useRemoveSampleTag(experimentId ?? 0, sampleId ?? 0);

  const exposureSummary = useMemo(
    () => ({
      total:    exposures.length,
      accepted: exposures.filter((e) => e.status === "accepted").length,
      rejected: exposures.filter((e) => e.status === "rejected").length,
    }),
    [exposures],
  );

  const handleSetStatus = useCallback(
    (status: "accepted" | "rejected" | null) => {
      if (!viewingId) return;
      setStatus.mutate({ exposureId: viewingId, status });
    },
    [viewingId, setStatus],
  );

  const handleSetIndexing = useCallback(() => {
    if (!viewingId) return;
    setIndexing.mutate(viewingId);
  }, [viewingId, setIndexing]);

  const handleAddTag = useCallback(
    (key: string, value: string) => {
      addExpTag.mutate({ key, value });
    },
    [addExpTag],
  );

  const handleRemoveTag = useCallback(
    (tagId: number) => {
      rmExpTag.mutate(tagId);
    },
    [rmExpTag],
  );

  if (!sample) return <div className="flex-1 min-h-0" />;

  return (
    <div
      data-testid="inspect-page"
      className="flex-1 min-h-0 flex flex-col gap-3 px-4 pb-4 pt-2"
    >
      <div className="flex items-center justify-end gap-2 px-2">
        <WarmAddMenu exposureId={viewingId} experimentId={experimentId} />
      </div>
      {/*
        Layout (shared with IndexPage via WorkspaceGrid):
          < 1400px  → single column stacked: image+gallery → metadata → chat
          ≥ 1400px  → three columns: chat | image+gallery | metadata
                      with minmax(320px,22fr) | 56fr | minmax(320px,22fr)
      */}
      <WorkspaceGrid
        left={<ChatCard />}
        center={
          <div className="flex flex-col gap-3 h-full min-h-0">
            <div className="flex-1 min-h-0">
              <Skeleton
                name="detector-image-card"
                className="h-full w-full"
                loading={exposuresQ.isLoading}
                stagger={50}
                transition={200}
                fixture={DETECTOR_IMAGE_FIXTURE}
                fallback={<div className="flex items-center justify-center h-full text-fg-muted text-sm">Loading exposure…</div>}
              >
                {viewingExposure ? (
                  <DetectorImageCard
                    exposure={viewingExposure}
                    onSetStatus={handleSetStatus}
                    onSetIndexing={handleSetIndexing}
                    onAddTag={handleAddTag}
                    onRemoveTag={handleRemoveTag}
                  />
                ) : (
                  <div className="flex items-center justify-center h-full text-fg-muted text-sm">
                    Select an exposure
                  </div>
                )}
              </Skeleton>
            </div>
            <div className="flex-none h-[140px] px-2 pt-3 pb-2 border-t border-border/40">
              <Skeleton
                name="thumbnail-gallery"
                className="h-full w-full"
                loading={exposuresQ.isLoading}
                stagger={50}
                transition={200}
                fixture={THUMBNAIL_GALLERY_FIXTURE}
                fallback={<div className="h-full w-full" />}
              >
                <ThumbnailGallery
                  exposures={exposures}
                  selectedId={viewingId}
                  onSelect={setViewingId}
                  className="h-full"
                />
              </Skeleton>
            </div>
          </div>
        }
        right={
          <Skeleton
            name="sample-metadata-card"
            className="h-full w-full"
            loading={samplesQ.isLoading}
            stagger={50}
            transition={200}
            fixture={SAMPLE_METADATA_FIXTURE}
            fallback={<div className="p-4 text-fg-muted text-sm">Loading sample…</div>}
          >
            <SampleMetadataCard
              sample={sample}
              experimentName={experimentName}
              exposureSummary={exposureSummary}
              onUpdateSample={(patch) => updateSample.mutate(patch)}
              onAddTag={(k, v) => addSampleTag.mutate({ key: k, value: v })}
              onRemoveTag={(id) => rmSampleTag.mutate(id)}
            />
          </Skeleton>
        }
        slotClassName={{
          left:   "min-h-[280px]",
          // Image-dominant card needs vertical room when stacked. The slot
          // is only this tall below 1400px; at the three-col breakpoint the
          // grid's fixed height takes over.
          center: "min-h-[640px]",
          right:  "min-h-[400px]",
        }}
      />
    </div>
  );
}
