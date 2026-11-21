import { Navigate, useParams } from "react-router-dom";
import { useExperiment } from "../../queries";
import { ApiError } from "../../api";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";

export function ExperimentEntry(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);
  if (exp.data === undefined && !exp.isError) return <ResolvingFallback />;
  if (exp.error instanceof ApiError && exp.error.status === 404) {
    return (
      <StaleUrlPage staleUrlContext={{
        kind: "not_found", missing: "experiment", missing_value: id ?? "",
        experiment_resolved: undefined, sample_resolved: undefined,
      }} />
    );
  }
  return <Navigate to={`/experiments/${expId}/workspace`} replace />;
}
