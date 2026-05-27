import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import { useExposures, usePeaks } from "../queries";
import { DetectorImage } from "./DetectorImage";
import { DetectorRingOverlay } from "./DetectorRingOverlay";
import { HintText } from "./ui";

/**
 * FocusDetectorPanel — the co-resident detector image on the focus workspace.
 *
 * Read-only image, with the q-link ring overlay layered on top (I4.3 / #180):
 * one ring per peak q, the ring matching `hoveredQ` lights, and hovering a
 * ring sets `hoveredQ` (cross-highlighting the trace peak). The ring overlay
 * is rotation-aware via the shared `decideOrient` rule.
 *
 * Reuses `DetectorImage` unchanged. We deliberately do NOT reuse the Inspect
 * `DetectorImageCard` — that bundles reject/tag/representative controls; the
 * focus surface's detector is a read-only companion to the trace.
 */
export function FocusDetectorPanel(): JSX.Element {
  const activeSampleId   = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const exposuresQ = useExposures(activeSampleId);
  const peaksQ = usePeaks(activeExposureId);
  const peakQs = (peaksQ.data ?? []).map((p) => p.q);

  const exposure = activeExposureId !== undefined
    ? exposuresQ.data?.find((e) => e.id === activeExposureId)
    : undefined;

  const body = (() => {
    if (activeSampleId === undefined) {
      return (
        <div data-testid="focus-detector-empty"
             className="flex-1 flex items-center justify-center">
          <HintText>Pick a sample to see its detector image.</HintText>
        </div>
      );
    }
    if (exposure === undefined || !exposure.image_path) {
      return (
        <div data-testid="focus-detector-empty"
             className="flex-1 flex items-center justify-center">
          <HintText>No detector image for this exposure.</HintText>
        </div>
      );
    }
    return (
      // Frame classes copied from LoupeFrame (the proven detector wrapper):
      // DetectorImage owns its own rotation via JS, so we only give it a box.
      <div className="relative mx-auto aspect-square w-full max-w-[420px]
                      overflow-hidden rounded border border-border bg-bg">
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          imageVersion={exposure.image_version}
          size="full"
          className="h-full w-full"
        />
        {peakQs.length > 0 && <DetectorRingOverlay peakQs={peakQs} />}
      </div>
    );
  })();

  return (
    // Card chrome mirrors the index cards: a bordered plate surface. `.panel`
    // is NOT an existing utility class — use explicit Tailwind/token classes.
    <section data-testid="focus-detector-panel"
             className="flex min-h-0 flex-col rounded border border-hair
                        bg-plate p-4">
      <div className="card-header">
        <span className="text-meta uppercase tracking-wider">Detector image</span>
      </div>
      <Skeleton
        name="focus-detector"
        className="flex-1 min-h-0 flex flex-col"
        loading={activeExposureId !== undefined && exposuresQ.isLoading}
        stagger={50}
        transition={200}
        fallback={
          <div className="flex-1 flex items-center justify-center">
            <HintText>Loading detector image…</HintText>
          </div>
        }
      >
        {body}
      </Skeleton>
    </section>
  );
}
