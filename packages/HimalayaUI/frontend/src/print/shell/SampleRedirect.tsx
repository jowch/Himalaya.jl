import { Navigate, useLocation, useParams } from "react-router-dom";
import { useCorpusSamples } from "../../queries";
import { ResolvingFallback } from "./ResolvingFallback";

export function workspaceHref(experimentId: number, sampleId: number, exposureId?: number | null): string {
  const base = `/experiments/${experimentId}/workspace/${sampleId}`;
  return exposureId ? `${base}?exposure=${exposureId}` : base;
}

export function SampleRedirect(): JSX.Element {
  const { sampleId } = useParams<{ sampleId: string }>();
  const { search } = useLocation();
  const corpus = useCorpusSamples();
  if (corpus.isLoading) return <ResolvingFallback />;
  const sample = (corpus.data ?? []).find((s) => s.id === Number(sampleId));
  const exposure = Number(new URLSearchParams(search).get("exposure"));
  return <Navigate to={sample ? workspaceHref(sample.experiment_id, sample.id, Number.isInteger(exposure) && exposure > 0 ? exposure : null) : "/experiments"} replace />;
}
