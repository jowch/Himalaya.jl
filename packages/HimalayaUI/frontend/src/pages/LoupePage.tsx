import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSamples,
  useExposures,
  usePeaks,
  useSetExposureStatus,
  useSelectExposure,
} from "../queries";
import { LoupeFrame } from "../components/LoupeFrame";
import { LoupeSidebar } from "../components/LoupeSidebar";
import type { Exposure } from "../api";

/**
 * Default exposure when the loupe opens: the indexing representative, else
 * the first accepted exposure, else the first exposure. Lifted from
 * InspectPage so the loupe opens on the same frame Inspect would have.
 */
function defaultExposureId(exposures: Exposure[]): number | undefined {
  const representative = exposures.find((e) => e.selected);
  if (representative) return representative.id;
  const firstAccepted = exposures.find((e) => e.status === "accepted");
  if (firstAccepted) return firstAccepted.id;
  return exposures[0]?.id;
}

// Fixture for boneyard's headless capture — a real render with mock props so
// the capture CLI can measure the loupe body layout (docs/boneyard.md Rule 2).
// image_path: null makes DetectorImage take its placeholder branch; the
// captured detector-frame bone is then a plain rectangle sized by LoupeFrame's
// aspect-square wrapper — correct geometry for a fixed-aspect frame.
const LOUPE_FIXTURE_EXPOSURE = {
  id: 0, sample_id: 0, filename: "JC000-001.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null,
  image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
const LOUPE_FIXTURE = (
  <div className="grid grid-cols-[1fr_286px] gap-7">
    <LoupeFrame
      exposure={LOUPE_FIXTURE_EXPOSURE}
      exposures={[LOUPE_FIXTURE_EXPOSURE]}
      onSelectExposure={() => {}}
    />
    <LoupeSidebar
      exposure={LOUPE_FIXTURE_EXPOSURE}
      exposures={[LOUPE_FIXTURE_EXPOSURE]}
      sample={{
        id: 0, experiment_id: 0, name: "JC000", display_name: "Sample",
        notes: null, tags: [], q_units: "A-1",
      }}
      signalLevel={0}
      onDropToggle={() => {}}
      onSetRepresentative={() => {}}
    />
  </div>
);

/**
 * LoupePage — the sample loupe at /samples/loupe/:sampleId. A focused
 * single-sample inspection surface: full detector image, exposure strip,
 * metadata sidebar. Mounted under the CorpusShell layout route (#161 / I1.5).
 *
 * URL-owned: the sample id comes from the route param, never from the
 * Zustand `activeSampleId` (master plan §2.3 — new surfaces own their URL).
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const navigate = useNavigate();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  const exposuresQ = useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);

  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);

  const isLoading = corpusQ.isLoading || exposuresQ.isLoading;

  // Active exposure — local state, defaulted by `defaultExposureId`. Reset
  // when the sample changes so the next sample picks its own default.
  const [activeId, setActiveId] = useState<number | undefined>(undefined);
  useEffect(() => {
    setActiveId(undefined);
  }, [sampleId]);
  const computedDefault = defaultExposureId(exposures);
  useEffect(() => {
    if (activeId === undefined && computedDefault !== undefined) {
      setActiveId(computedDefault);
    }
  }, [activeId, computedDefault]);
  const activeExposure = exposures.find((e) => e.id === activeId);

  // M-8: signal-strength meter level — an honest proxy from the active
  // exposure's detected-peak count (more resolvable Bragg reflections ⇒
  // stronger ordering signal). The meter clamps to 0..5.
  const peaksQ = usePeaks(activeExposure?.id);
  const signalLevel = peaksQ.data?.length ?? 0;

  // L-9: subtitle reads `<sample-id> · exposure N of M`. Position is the
  // active exposure's 1-based index in the per-sample exposure list.
  const frameIndex = exposures.findIndex((e) => e.id === activeId);
  const exposurePosition =
    frameIndex >= 0 ? `exposure ${frameIndex + 1} of ${exposures.length}` : "—";

  const setStatus = useSetExposureStatus(hasValidId ? sampleId : 0);
  const setRepresentative = useSelectExposure(hasValidId ? sampleId : 0);

  const handleDropToggle = useCallback(() => {
    if (!activeExposure) return;
    setStatus.mutate({
      exposureId: activeExposure.id,
      status: activeExposure.status === "rejected" ? null : "rejected",
    });
  }, [activeExposure, setStatus]);

  const handleSetRepresentative = useCallback(() => {
    if (!activeExposure) return;
    setRepresentative.mutate(activeExposure.id);
  }, [activeExposure, setRepresentative]);

  const flip = useCallback(
    (delta: number) => {
      if (activeId === undefined || exposures.length === 0) return;
      const idx = exposures.findIndex((e) => e.id === activeId);
      if (idx < 0) return;
      const next = Math.min(Math.max(idx + delta, 0), exposures.length - 1);
      setActiveId(exposures[next].id);
    },
    [activeId, exposures],
  );

  const goBack = useCallback(() => {
    navigate("/samples");
  }, [navigate]);

  // Loupe keyboard shortcuts. The input/textarea guard is in from the start
  // so it survives #159 adding tag-edit inputs to the sidebar.
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      const tag = (e.target as HTMLElement | null)?.tagName;
      if (tag === "INPUT" || tag === "TEXTAREA") return;
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
      <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
        <div
          data-testid="loupe-not-found"
          className="rounded border border-border p-8 text-sm text-ink-faint"
        >
          Sample not found.{" "}
          <button
            onClick={goBack}
            className="font-semibold text-accent hover:underline"
          >
            Back to the sheet
          </button>
        </div>
      </div>
    );
  }

  return (
    <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
      <button
        data-testid="loupe-back"
        onClick={goBack}
        className="mb-3.5 text-[11.5px] font-semibold text-accent hover:underline"
      >
        ← Back to the sheet
      </button>
      <header className="mb-5">
        <h2 className="text-2xl text-ink">
          {sample?.display_name ?? sample?.name ?? "—"}
        </h2>
        <p
          data-testid="loupe-subtitle"
          className="mt-1 font-mono text-xs text-ink-faint"
        >
          {sample?.name ?? "—"} · {exposurePosition}
        </p>
      </header>

      <Skeleton
        name="loupe"
        className="block"
        loading={isLoading}
        stagger={50}
        transition={200}
        fixture={LOUPE_FIXTURE}
        fallback={
          <div
            data-testid="loupe-skeleton"
            className="p-8 text-sm italic text-ink-faint"
          >
            Loading sample…
          </div>
        }
      >
        <div className="grid grid-cols-[1fr_286px] gap-7">
          {sample && activeExposure ? (
            <>
              <LoupeFrame
                exposure={activeExposure}
                exposures={exposures}
                onSelectExposure={setActiveId}
              />
              <LoupeSidebar
                exposure={activeExposure}
                exposures={exposures}
                sample={sample}
                signalLevel={signalLevel}
                onDropToggle={handleDropToggle}
                onSetRepresentative={handleSetRepresentative}
              />
            </>
          ) : (
            <div className="col-span-2 p-8 text-sm text-ink-faint">
              This sample has no exposures.
            </div>
          )}
        </div>
      </Skeleton>
    </div>
  );
}
