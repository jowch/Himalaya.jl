import { useEffect, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { useQueryClient, useQuery } from "@tanstack/react-query";
import { useAppState } from "../state";
import { useExperiments, queryKeys } from "../queries";
import * as api from "../api";
import type { Experiment, Sample, Exposure } from "../api";
import { consumeEmitMode } from "../lib/url/emitMode";

// Spec §4.3 — Zustand → URL via react-router useNavigate. Subscribes to
// experiments + samples queries via TanStack so SSE-driven cache rewrites
// trigger a re-render of this hook (§7 invalidation).

function nameForExperiment(list: Experiment[] | undefined, id: number | undefined):
  string | undefined
{
  if (id === undefined || list === undefined) return undefined;
  const found = list.find((e) => e.id === id);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function nameForSample(list: Sample[] | undefined, sId: number | undefined):
  string | undefined
{
  if (sId === undefined || list === undefined) return undefined;
  const found = list.find((s) => s.id === sId);
  if (found === undefined) return undefined;
  return found.name === null ? undefined : found.name;
}

function filenameForExposure(qc: ReturnType<typeof useQueryClient>,
                              sId: number | undefined,
                              eId: number | undefined): string | undefined
{
  if (sId === undefined || eId === undefined) return undefined;
  const list = qc.getQueryData<Exposure[]>(queryKeys.exposures(sId)) ?? [];
  const found = list.find((e) => e.id === eId);
  if (found === undefined) return undefined;
  return found.filename === null ? undefined : found.filename;
}

function buildUrl(
  page: "index" | "inspect" | "compare",
  experiment: string | undefined,
  sample: string | undefined,
  exposure: string | undefined,
  current: string,
): string {
  const enc = (s: string) => encodeURIComponent(s);
  if (page === "compare") {
    // Don't try to re-emit a Compare URL — Compare uses numeric ids that
    // useStateFromUrl doesn't track in Zustand. ComparePage handles its
    // own URL via useNavigate. Returning current keeps the equality guard
    // happy and prevents accidental redirect.
    return current;
  }
  const parts = [`/${page}`];
  if (experiment !== undefined) {
    parts.push(enc(experiment));
    if (sample !== undefined) parts.push(enc(sample));
  }
  let url = parts.join("/");
  if (page === "inspect" && exposure !== undefined && sample !== undefined) {
    url += `?exposure=${enc(exposure)}`;
  }
  return url;
}

export function useUrlFromState(): void {
  const navigate = useNavigate();
  const location = useLocation();
  const qc = useQueryClient();

  const activePage = useAppState((s) => s.activePage);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);

  // Subscribe to experiments + (when an experiment is active) samples
  // queries so SSE-driven cache rewrites trigger a re-render. We use
  // `useQuery` directly here (rather than `useSamples`) so we can gate
  // on `enabled: activeExperimentId !== undefined` — `useSamples` doesn't
  // expose `enabled`, and calling it with `0` would fire `GET
  // /api/experiments/0/samples` → 404 with retries on every cold mount.
  const { data: experiments } = useExperiments();
  const samplesQuery = useQuery({
    queryKey: queryKeys.samples(activeExperimentId ?? 0),
    queryFn: () => api.listSamples(activeExperimentId as number),
    enabled: activeExperimentId !== undefined,
  });
  const samples = activeExperimentId !== undefined ? samplesQuery.data : undefined;

  // Track the previous resolved slug pair so we can detect SSE-driven
  // disappearance (a slug went from defined → undefined because the
  // entity was deleted from the cache). Per spec §4.3 + §7, that case
  // should emit replace, not push (otherwise back-button stops at the
  // broken URL).
  const prevSlugsRef = useRef<{ exp: string | undefined; sample: string | undefined }>(
    { exp: undefined, sample: undefined },
  );

  useEffect(() => {
    const expName = nameForExperiment(experiments, activeExperimentId);
    const sampleName = nameForSample(samples, activeSampleId);
    const exposureName = filenameForExposure(qc, activeSampleId, activeExposureId);
    const current = location.pathname + location.search;
    const target = buildUrl(activePage, expName, sampleName, exposureName, current);

    // SSE-driven invalidation detection: a previously-resolvable slug
    // is now undefined (the row vanished from cache). Force replace.
    const prev = prevSlugsRef.current;
    const slugDisappeared =
      (prev.exp !== undefined && expName === undefined && activeExperimentId !== undefined) ||
      (prev.sample !== undefined && sampleName === undefined && activeSampleId !== undefined);
    prevSlugsRef.current = { exp: expName, sample: sampleName };

    if (target === current) return;        // equality guard
    const explicitMode = consumeEmitMode();
    const mode = slugDisappeared ? "replace" : explicitMode;
    navigate(target, { replace: mode === "replace" });
  }, [
    activePage, activeExperimentId, activeSampleId, activeExposureId,
    experiments, samples,
    location.pathname, location.search,
    navigate, qc,
  ]);
}
