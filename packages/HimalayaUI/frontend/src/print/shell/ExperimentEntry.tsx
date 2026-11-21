import { Navigate, useParams } from "react-router-dom";
import { useExperiment } from "../../queries";
import { ResolvingFallback } from "./ResolvingFallback";

export function ExperimentEntry(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);
  if (exp.data === undefined && !exp.isError) return <ResolvingFallback />;
  return <Navigate to={`/experiments/${expId}/workspace`} replace />;
}
