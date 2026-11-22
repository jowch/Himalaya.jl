import { Navigate, useParams } from "react-router-dom";
import { useExperiment } from "../../queries";
import { ApiError } from "../../api";
import { ResolvingFallback } from "./ResolvingFallback";
import { LinkGone } from "./LinkGone";

export const EXPERIMENT_GONE = "That experiment no longer exists.";

export function ExperimentEntry(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);
  if (exp.data === undefined && !exp.isError) return <ResolvingFallback />;
  if (exp.error instanceof ApiError && exp.error.status === 404) return <LinkGone message={EXPERIMENT_GONE} />;
  return <Navigate to={`/experiments/${expId}/workspace`} replace />;
}
