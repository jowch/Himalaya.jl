import { useParams } from "react-router-dom";
import { useCorpusSamples, useExposures, useExperiment } from "../queries";

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
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const experimentQ = useExperiment(sample?.experiment_id ?? 0);
  const experimentName =
    experimentQ.data?.name ?? experimentQ.data?.path ?? undefined;

  return (
    <div
      data-testid="loupe-page"
      className="mx-auto max-w-[1100px] px-8 py-7"
    >
      <header className="mb-5">
        <h2 className="text-2xl text-ink">
          {sample?.display_name ?? sample?.name ?? "—"}
        </h2>
        <p className="mt-1 font-mono text-xs text-ink-faint">
          {experimentName ?? "—"}
          {sample?.name ? ` · ${sample.name}` : ""}
        </p>
      </header>
    </div>
  );
}
