import { useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { useCorpusSamples, useExposures, useExperiment } from "../queries";
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
  const experimentQ = useExperiment(sample?.experiment_id ?? 0);
  const experimentName =
    experimentQ.data?.name ?? experimentQ.data?.path ?? undefined;

  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);

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

  function goBack(): void {
    navigate("/samples");
  }

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
        <p className="mt-1 font-mono text-xs text-ink-faint">
          {experimentName ?? "—"}
          {sample?.name ? ` · ${sample.name}` : ""}
        </p>
      </header>

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
              onDropToggle={() => {}}
              onSetRepresentative={() => {}}
            />
          </>
        ) : (
          <div className="col-span-2 p-8 text-sm text-ink-faint">
            This sample has no exposures.
          </div>
        )}
      </div>
    </div>
  );
}
