import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import { useExposures, usePeaks, useSelectExposure } from "../queries";
import { DetectorImage } from "./DetectorImage";
import { DetectorRingOverlay } from "./DetectorRingOverlay";
import { Card, HintText } from "./ui";

/**
 * FocusDetectorPanel — the co-resident detector image on the focus workspace.
 *
 * Read-only image, with the q-link ring overlay layered on top (I4.3 / #180):
 * one ring per peak q, the ring matching `hoveredQ` lights, and hovering a
 * ring sets `hoveredQ` (cross-highlighting the trace peak). The ring overlay
 * is rotation-aware via the shared `decideOrient` rule.
 *
 * Reuses `DetectorImage` unchanged. We deliberately do NOT reuse the Inspect
 * `DetectorImageCard` — that bundles reject/tag controls; the focus surface's
 * detector is a read-only companion to the trace.
 *
 * R5 (#228, F-11): the panel carries a representative-exposure switcher — a
 * strip of mini DetectorImage thumbnails, one per exposure. Clicking a thumb
 * switches the VIEWED exposure (local view state via `setActiveExposure`);
 * "Set representative" PERSISTS the choice (`useSelectExposure` →
 * `PATCH /exposures/:id/select`, the same mechanism the Inspect loupe uses).
 * The persisted representative is the `selected` exposure (rep marker). The
 * strip is suppressed when a sample has only one exposure (nothing to switch).
 */
export function FocusDetectorPanel(): JSX.Element {
  const activeSampleId    = useAppState((s) => s.activeSampleId);
  const activeExposureId  = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  const exposuresQ = useExposures(activeSampleId);
  const setRepresentative = useSelectExposure(activeSampleId ?? 0);
  // Deliberately NOT gated by `useExposureHasPendingPeakOps` (src/AGENTS.md):
  // the rings are presentational and keyed by q, so a transient optimistic
  // peak landing mid-mutation just draws a brief extra ring — there's no
  // derived/aggregate value to hide. Gating here would needlessly blank the
  // q-link rings during every peak edit. Left ungated on purpose.
  const peaksQ = usePeaks(activeExposureId);
  const peakQs = (peaksQ.data ?? []).map((p) => p.q);

  const exposures = exposuresQ.data ?? [];
  const exposure = activeExposureId !== undefined
    ? exposures.find((e) => e.id === activeExposureId)
    : undefined;

  // F-11: the switcher only earns its space when there's more than one
  // exposure to switch between. The viewed exposure is `exposure`; the
  // persisted representative is whichever exposure has `selected: true`.
  const showSwitcher = exposures.length > 1;
  const viewedIsRep = exposure?.selected === true;

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
      // T-8 (R0c): the detector frame is the dark `frame-edge` window, not paper.
      <div className="relative mx-auto aspect-square w-full max-w-[420px]
                      overflow-hidden rounded border border-frame-edge bg-frame-edge">
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
    // R3-N3 (#209): the detector panel is plate-like and belongs to the
    // "Plate Lift" family — DESIGN.md §Elevation says the figure plate and
    // plate-like cards are the one elevated object. Promote the explicit
    // `rounded border bg-plate` triple to the canonical `.card` rule so the
    // panel floats above the paper with the same warm shadow as the trace
    // plate. Keep `p-4` for body padding.
    <Card as="section" elevated data-testid="focus-detector-panel"
             className="flex min-h-0 flex-col p-4">
      <div className="card-header flex items-center justify-between gap-3">
        <span className="text-meta uppercase tracking-wider">Detector image</span>
        {showSwitcher && (
          <div
            data-testid="exposure-switcher"
            className="flex items-center gap-1.5"
          >
            {exposures.map((e) => {
              const isViewed = e.id === activeExposureId;
              return (
                <button
                  key={e.id}
                  type="button"
                  data-testid={`exposure-thumb-${e.id}`}
                  data-rep={e.selected ? "true" : undefined}
                  data-viewed={isViewed ? "true" : undefined}
                  title={
                    (e.filename ?? `exposure ${e.id}`) +
                    (e.selected ? " (representative)" : "")
                  }
                  onClick={() => setActiveExposure(e.id)}
                  className={
                    "relative h-8 w-8 overflow-hidden rounded-sm border " +
                    "border-frame-edge bg-frame-edge " +
                    (isViewed ? "ring-2 ring-print-accent" : "")
                  }
                >
                  {e.image_path ? (
                    <DetectorImage
                      exposureId={e.id}
                      imagePath={e.image_path}
                      imageVersion={e.image_version}
                      size="thumb"
                      className="h-full w-full"
                    />
                  ) : (
                    <span className="block h-full w-full" />
                  )}
                  {e.selected && (
                    <span
                      aria-hidden="true"
                      className="absolute right-0.5 top-0.5 h-1.5 w-1.5
                                 rounded-full border border-plate bg-print-accent"
                    />
                  )}
                </button>
              );
            })}
            <button
              type="button"
              data-testid="exposure-set-rep"
              disabled={viewedIsRep || exposure === undefined}
              onClick={() => exposure && setRepresentative.mutate(exposure.id)}
              className="rounded border border-hair-strong px-1.5 py-0.5 text-[10px]
                         font-semibold uppercase tracking-wide text-ink
                         disabled:cursor-not-allowed disabled:opacity-40"
            >
              {viewedIsRep ? "Representative" : "Set rep"}
            </button>
          </div>
        )}
      </div>
      <Skeleton
        name="focus-detector"
        className="flex-1 min-h-0 flex flex-col"
        loading={activeExposureId !== undefined && exposuresQ.isLoading}
        // U-2 (#255): a detector window is dark even while loading, so its
        // skeleton shimmers `frame-edge` (the dark window backing), overriding
        // the global `paper-sunk` bone fill used on light surfaces. oklch literal
        // matching `--color-frame-edge` (boneyard takes a plain CSS color string).
        color="oklch(0.150 0.010 55)"
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
    </Card>
  );
}
