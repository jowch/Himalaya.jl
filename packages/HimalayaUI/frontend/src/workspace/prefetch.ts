import { useEffect } from "react";
import { useQueryClient, type QueryClient } from "@tanstack/react-query";
import * as api from "../api";
import type { CorpusSample, Exposure } from "../api";
import { queryKeys } from "../queries";
import { imageUrl, nextUnreviewed } from "./model";

const preloading = new Set<HTMLImageElement>();

export function preloadImage(url: string, priority: "high" | "low" = "low"): void {
  if (typeof Image === "undefined") return;
  const img = new Image();
  img.fetchPriority = priority;
  img.decoding = "async";
  const done = (): void => {
    preloading.delete(img);
  };
  img.onload = done;
  img.onerror = done;
  preloading.add(img);
  img.src = url;
}

export function prefetchExposureData(qc: QueryClient, exposureId: number): Promise<unknown> {
  return Promise.all([
    qc.prefetchQuery({ queryKey: queryKeys.trace(exposureId), queryFn: () => api.getTrace(exposureId) }),
    qc.prefetchQuery({ queryKey: queryKeys.peaks(exposureId), queryFn: () => api.listPeaks(exposureId) }),
    qc.prefetchQuery({ queryKey: queryKeys.indices(exposureId), queryFn: () => api.listIndices(exposureId) }),
    qc.prefetchQuery({ queryKey: queryKeys.assignment(exposureId), queryFn: () => api.getAssignment(exposureId) }),
  ]);
}

export async function prefetchSample(
  qc: QueryClient,
  sample: Pick<CorpusSample, "id" | "representative_exposure_id">,
  { fullImages }: { fullImages: boolean },
): Promise<void> {
  const rep = sample.representative_exposure_id ?? undefined;
  const data = rep === undefined ? Promise.resolve() : prefetchExposureData(qc, rep);
  await qc.prefetchQuery({ queryKey: queryKeys.exposures(sample.id), queryFn: () => api.listExposures(sample.id) });
  const exposures = qc.getQueryData<Exposure[]>(queryKeys.exposures(sample.id)) ?? [];
  for (const e of exposures) {
    const thumb = imageUrl(e, true);
    if (thumb) preloadImage(thumb);
  }
  await data;
  if (!fullImages) return;
  for (const e of exposures) {
    const full = imageUrl(e, false);
    if (full) preloadImage(full);
  }
}

export function neighbourSamples(samples: readonly CorpusSample[], currentId: number | undefined): CorpusSample[] {
  const at = samples.findIndex((s) => s.id === currentId);
  if (at < 0) return [];
  const ids = new Set<number>();
  const prev = samples[at - 1];
  const next = samples[at + 1];
  if (next) ids.add(next.id);
  if (prev) ids.add(prev.id);
  const unreviewed = nextUnreviewed(samples, currentId, 1);
  if (unreviewed !== undefined) ids.add(unreviewed);
  ids.delete(samples[at]!.id);
  return samples.filter((s) => ids.has(s.id));
}

export function usePrefetchNeighbours(
  samples: readonly CorpusSample[],
  current: CorpusSample | undefined,
  currentReady: boolean,
): void {
  const qc = useQueryClient();
  useEffect(() => {
    if (!current || !currentReady) return;
    void prefetchSample(qc, current, { fullImages: true });
    for (const s of neighbourSamples(samples, current.id)) void prefetchSample(qc, s, { fullImages: false });
  }, [qc, samples, current, currentReady]);
}
