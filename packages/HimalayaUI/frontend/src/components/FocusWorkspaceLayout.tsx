import { useAppState } from "../state";
import {
  useCorpusSamples, useSamples, useUpdateSample,
  useExperiment, useExposures,
} from "../queries";
import { sampleDisplayName } from "../lib/sample/displayName";
import { PlotCard } from "./PlotCard";
import { FocusPlotHeader } from "./FocusPlotHeader";
import { IndicesCard } from "./IndicesCard";
import { FocusDetectorPanel } from "./FocusDetectorPanel";
import { FocusNotesMargin } from "./FocusNotesMargin";

/**
 * FocusWorkspaceLayout — the focus-workspace grid (mockup focus-workspace.html):
 *
 *   [ work (trace hero + co-resident detector) | rail (phase call) | notes ]
 *
 * Pure composition. Every panel is an existing component reused unchanged:
 *   - trace hero  = <PlotCard/>      (also drives useAutoPickExposure)
 *   - rail        = <IndicesCard/>   (PhasePanel phase call + candidates + Miller)
 *   - detector    = <FocusDetectorPanel/> (read-only DetectorImage; q-link is I4.3)
 *   - notes       = <FocusNotesMargin/>
 *
 * The notes margin is the first column to yield: it collapses below the
 * `xl` breakpoint (mirrors the mockup's max-width:1320px rule). The q-link
 * cross-highlight (#180) and the Index cutover (#181) are out of scope.
 */
export function FocusWorkspaceLayout(): JSX.Element {
  const activeSampleId   = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);

  // The corpus query gives us the sample's authoritative `experiment_id`
  // (the I4.1 route shim seeds only `activeSampleId`, NOT activeExperimentId).
  const corpusQ = useCorpusSamples();
  const corpusSample = activeSampleId !== undefined
    ? corpusQ.data?.find((s) => s.id === activeSampleId)
    : undefined;
  const experimentId = corpusSample?.experiment_id;

  // ── Focus plate header (R3 / #226) ────────────────────────────────────────
  // The trace plate's header is seeded from the ROUTE's sample (not the global
  // experiment-picker, which isn't seeded on /sample/:sampleId). We resolve:
  //   serif title  = display_name || name || "Sample N"  (corpusSample)
  //   sub: code     = sample.name (the internal code, e.g. "smp_09")
  //   sub: beamtime = owning experiment name
  //   sub: exposure = the active (or selected) exposure's filename stem
  const experimentQ = useExperiment(experimentId ?? 0);
  const exposuresQ  = useExposures(activeSampleId);
  const repExposure = activeExposureId !== undefined
    ? exposuresQ.data?.find((e) => e.id === activeExposureId)
    : exposuresQ.data?.find((e) => e.selected);
  const exposureLabel = repExposure?.filename
    ? repExposure.filename.replace(/\.[^.]+$/, "")
    : null;
  const focusHeader = corpusSample !== undefined
    ? (
      <FocusPlotHeader
        sampleName={sampleDisplayName(corpusSample)}
        sampleCode={corpusSample.name}
        beamtime={experimentQ.data?.name ?? null}
        exposureLabel={exposureLabel}
      />
    )
    : undefined;

  // CACHE-COHERENCE (blocking-review fix): the notes textarea must read from
  // the SAME cache `updateSampleMutator` patches. That mutator writes
  // queryKeys.samples(experimentId) + queryKeys.sample(sampleId) — it does
  // NOT touch queryKeys.corpusSamples. So we source the notes sample from
  // `useSamples(experimentId)` (the experiment-scoped list), exactly as the
  // proven Inspect path does (InspectPage.tsx:125,161 -> SampleMetadataCard ->
  // useUpdateSample). Reading notes off `useCorpusSamples` instead would make
  // the on-blur focus-gate re-sync a stale value and silently revert the edit.
  const samplesQ = useSamples(experimentId ?? 0);
  const notesSample = activeSampleId !== undefined && experimentId !== undefined
    ? samplesQ.data?.find((s) => s.id === activeSampleId)
    : undefined;

  // Write target uses the authoritative experiment_id (NOT the nullable
  // activeExperimentId, which isn't seeded on /sample/:sampleId). Hook is
  // always called (hooks rule); it only fires from the margin's blur.
  const updateSample = useUpdateSample(experimentId ?? 0, activeSampleId ?? 0);

  return (
    <div
      data-testid="focus-workspace-layout"
      className="grid min-h-0 flex-1 grid-cols-1 xl:grid-cols-[1fr_348px_250px]"
    >
      {/* work area: trace hero stacked over the co-resident detector */}
      <div className="flex min-h-0 flex-col gap-5 overflow-auto p-6">
        <div className="min-h-[420px]">
          <PlotCard {...(focusHeader ? { headerSlot: focusHeader } : {})} />
        </div>
        <FocusDetectorPanel />
      </div>

      {/* rail: the phase call + candidates */}
      <aside className="min-h-0 overflow-y-auto border-l border-hair bg-paper-sunk">
        <IndicesCard />
      </aside>

      {/* notes margin: collapses below xl. Gated on the experiment-scoped
          sample being loaded so the textarea never reads a stale source. */}
      {notesSample !== undefined && (
        <div className="hidden xl:block">
          <FocusNotesMargin
            sample={notesSample}
            onSaveNotes={(notes) => updateSample.mutate({ notes })}
          />
        </div>
      )}
    </div>
  );
}
