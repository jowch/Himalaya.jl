import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useAppState } from "../state";
import {
  useExperiment,
  useExposures,
  useSamples,
  useSetExposureStatus,
  useSelectExposure,
  useAddExposureTag,
  useAddSampleTag,
  useRemoveSampleTag,
  useUpdateSample,
} from "../queries";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { DetectorImageCard } from "../components/DetectorImageCard";
import { SampleMetadataCard } from "../components/SampleMetadataCard";
import { ChatCard } from "../components/ChatCard";
import { WorkspaceGrid } from "../components/WorkspaceGrid";

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
  useEffect(() => {
    if (viewingId === undefined && defaultId !== undefined) {
      setViewingId(defaultId);
    }
  }, [defaultId, viewingId]);

  const viewingExposure = exposures.find((e) => e.id === viewingId);

  const setStatus    = useSetExposureStatus(sampleId ?? 0);
  const setIndexing  = useSelectExposure(sampleId ?? 0);
  const addExpTag    = useAddExposureTag(sampleId ?? 0, viewingId ?? 0);
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

  if (!sample) return <div className="flex-1 min-h-0" />;

  return (
    <div
      data-testid="inspect-page"
      className="flex-1 min-h-0 flex flex-col gap-3 px-4 pb-6 pt-2"
    >
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
              {viewingExposure ? (
                <DetectorImageCard
                  exposure={viewingExposure}
                  onSetStatus={handleSetStatus}
                  onSetIndexing={handleSetIndexing}
                  onAddTag={handleAddTag}
                />
              ) : (
                <div className="flex items-center justify-center h-full text-fg-muted text-sm">
                  Select an exposure
                </div>
              )}
            </div>
            <div className="flex-none h-[140px] px-2 pt-3 pb-2 border-t border-border/40">
              <ThumbnailGallery
                exposures={exposures}
                selectedId={viewingId}
                onSelect={setViewingId}
                className="h-full"
              />
            </div>
          </div>
        }
        right={
          <SampleMetadataCard
            sample={sample}
            experimentName={experimentName}
            exposureSummary={exposureSummary}
            onUpdateSample={(patch) => updateSample.mutate(patch)}
            onAddTag={(k, v) => addSampleTag.mutate({ key: k, value: v })}
            onRemoveTag={(id) => rmSampleTag.mutate(id)}
          />
        }
        slotClassName={{
          left:   "min-h-[280px]",
          // Image-dominant card needs vertical room when stacked. The slot
          // is only this tall below 1400px; at the three-col breakpoint the
          // grid's fixed height takes over.
          center: "min-h-[640px]",
          right:  "min-h-[200px]",
        }}
      />
    </div>
  );
}
