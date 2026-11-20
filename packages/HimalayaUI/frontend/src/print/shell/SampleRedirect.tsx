import { Navigate, useParams } from "react-router-dom";
import { useCorpusSamples } from "../../queries";
import { ResolvingFallback } from "./ResolvingFallback";

export function SampleRedirect(): JSX.Element {
  const { sampleId } = useParams<{ sampleId: string }>();
  const corpus = useCorpusSamples();
  if (corpus.isLoading) return <ResolvingFallback />;
  const sample = (corpus.data ?? []).find((s) => s.id === Number(sampleId));
  return <Navigate to={sample ? `/experiments/${sample.experiment_id}/workspace/${sample.id}` : "/experiments"} replace />;
}
