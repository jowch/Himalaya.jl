import { useEffect } from "react";
import { Navigate, useLocation, useParams } from "react-router-dom";
import { showToast } from "../../lib/toast";
import { useCorpusSamples } from "../../queries";
import { ResolvingFallback } from "./ResolvingFallback";

export const STALE_LINK = "That sample link is no longer valid";

export function workspaceHref(experimentId: number, sampleId: number, exposureId?: number | null): string {
  const base = `/experiments/${experimentId}/workspace/${sampleId}`;
  return exposureId ? `${base}?exposure=${exposureId}` : base;
}

export function SampleRedirect(): JSX.Element {
  const { sampleId } = useParams<{ sampleId: string }>();
  const { search } = useLocation();
  const corpus = useCorpusSamples();
  const sample = (corpus.data ?? []).find((s) => s.id === Number(sampleId));
  const stale = !corpus.isLoading && !sample;
  useEffect(() => {
    if (stale) showToast(STALE_LINK, "warning");
  }, [stale]);
  if (corpus.isLoading) return <ResolvingFallback />;
  const exposure = Number(new URLSearchParams(search).get("exposure"));
  return <Navigate to={sample ? workspaceHref(sample.experiment_id, sample.id, Number.isInteger(exposure) && exposure > 0 ? exposure : null) : "/experiments"} replace />;
}
