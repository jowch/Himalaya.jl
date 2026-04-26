import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useAppState } from "../state";
import {
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

  const exposuresQ = useExposures(sampleId);
  const samplesQ   = useSamples(experimentId ?? 0);
  const exposures  = exposuresQ.data ?? [];
  const sample     = samplesQ.data?.find((s) => s.id === sampleId);

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
      className="flex-1 min-h-0 flex flex-col gap-3 px-4 pb-4 pt-2"
    >
      {/*
        Breakpoints:
          < 1100px  → single column stacked
          1100–1400 → two columns 3fr/2fr: left=(metadata+gallery) right=image
          ≥ 1400px  → three columns 28fr/22fr/50fr: metadata | gallery | image
      */}
      <div
        className="
          min-h-0 grid gap-3 flex-1
          grid-cols-1
          min-[1100px]:grid-cols-[3fr_2fr]
          min-[1400px]:grid-cols-[28fr_22fr_50fr]
          h-auto min-[1100px]:h-[700px]
          min-[1100px]:max-h-[calc(100vh-120px)]
        "
      >
        {/* Metadata — col 1 row 1 on medium, col 1 on large */}
        <section
          className="
            card min-h-[200px] min-[1100px]:min-h-0 overflow-hidden
            order-1
            min-[1100px]:row-start-1 min-[1100px]:col-start-1
            min-[1400px]:col-start-1
          "
        >
          <SampleMetadataCard
            sample={sample}
            exposureSummary={exposureSummary}
            onUpdateSample={(patch) => updateSample.mutate(patch)}
            onAddTag={(k, v) => addSampleTag.mutate({ key: k, value: v })}
            onRemoveTag={(id) => rmSampleTag.mutate(id)}
          />
        </section>

        {/* Thumbnail gallery — col 1 row 2 on medium, col 2 on large */}
        <section
          className="
            card min-h-[180px] min-[1100px]:min-h-0 overflow-hidden p-2
            order-2
            min-[1100px]:row-start-2 min-[1100px]:col-start-1
            min-[1400px]:row-start-1 min-[1400px]:col-start-2
          "
        >
          <ThumbnailGallery
            exposures={exposures}
            selectedId={viewingId}
            onSelect={setViewingId}
            columns="auto"
            className="h-full"
          />
        </section>

        {/* Detector image — col 2 row-span-2 on medium, col 3 on large */}
        <section
          className="
            card min-h-[300px] min-[1100px]:min-h-0 overflow-hidden
            order-3
            min-[1100px]:row-span-2 min-[1100px]:col-start-2
            min-[1400px]:row-span-1 min-[1400px]:col-start-3
          "
        >
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
        </section>
      </div>
    </div>
  );
}
