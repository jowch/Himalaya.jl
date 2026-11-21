import { useEffect, useRef } from "react";
import { useQueryClient, type QueryClient, type QueryKey } from "@tanstack/react-query";
import * as api from "../api";
import type { CorpusSample, Exposure } from "../api";
import { queryKeys } from "../queries";
import { imageUrl, nextUnreviewed } from "./model";

export const PREFETCH_DWELL_MS = 150;

export function preloadImage(url: string, signal: AbortSignal): Promise<void> {
  if (typeof Image === "undefined" || signal.aborted) return Promise.resolve();
  return new Promise((resolve) => {
    const img = new Image();
    img.decoding = "async";
    const finish = (): void => {
      signal.removeEventListener("abort", abort);
      resolve();
    };
    const abort = (): void => {
      img.onload = null;
      img.onerror = null;
      img.src = "";
      finish();
    };
    img.onload = finish;
    img.onerror = finish;
    signal.addEventListener("abort", abort);
    img.src = url;
  });
}

type Neighbour = { id: number; rep: number | undefined };

class PrefetchRun {
  readonly controller = new AbortController();
  readonly keys: QueryKey[] = [];

  constructor(private readonly qc: QueryClient) {}

  get signal(): AbortSignal {
    return this.controller.signal;
  }

  query<T>(queryKey: QueryKey, fetcher: (signal: AbortSignal) => Promise<T>): Promise<void> {
    if (this.signal.aborted) return Promise.resolve();
    this.keys.push(queryKey);
    return this.qc.prefetchQuery({ queryKey, queryFn: ({ signal }) => fetcher(signal) });
  }

  exposureData(id: number): Promise<unknown> {
    return Promise.all([
      this.query(queryKeys.trace(id), (signal) => api.getTrace(id, { signal })),
      this.query(queryKeys.peaks(id), (signal) => api.listPeaks(id, { signal })),
      this.query(queryKeys.indices(id), (signal) => api.listIndices(id, { signal })),
      this.query(queryKeys.assignment(id), (signal) => api.getAssignment(id, { signal })),
    ]);
  }

  async neighbour({ id, rep }: Neighbour): Promise<void> {
    const data = rep === undefined ? Promise.resolve() : this.exposureData(rep);
    await this.query(queryKeys.exposures(id), (signal) => api.listExposures(id, { signal }));
    const thumbs = this.exposures(id).flatMap((e) => imageUrl(e, true) ?? []);
    await Promise.all([data, ...thumbs.map((url) => preloadImage(url, this.signal))]);
  }

  exposures(sampleId: number): Exposure[] {
    return this.qc.getQueryData<Exposure[]>(queryKeys.exposures(sampleId)) ?? [];
  }

  async images(urls: readonly string[]): Promise<void> {
    for (const url of urls) {
      if (this.signal.aborted) return;
      await preloadImage(url, this.signal);
    }
  }

  cancel(): void {
    this.controller.abort();
    const keys = this.keys.splice(0);
    setTimeout(() => {
      for (const queryKey of keys) {
        const query = this.qc.getQueryCache().find({ queryKey, exact: true });
        if (query && query.state.fetchStatus === "fetching" && query.getObserversCount() === 0) {
          void this.qc.cancelQueries({ queryKey, exact: true });
        }
      }
    }, 0);
  }
}

async function prefetchAround(run: PrefetchRun, shownFull: string[], otherFull: string[], neighbours: Neighbour[]): Promise<void> {
  await Promise.all([run.images(shownFull), ...neighbours.map((n) => run.neighbour(n))]);
  await run.images(otherFull);
  await run.images(
    neighbours.flatMap((n) => {
      const rep = run.exposures(n.id).find((e) => e.id === n.rep);
      return (rep && imageUrl(rep, false)) ?? [];
    }),
  );
}

export function neighbourSamples(samples: readonly CorpusSample[], currentId: number | undefined): CorpusSample[] {
  const at = samples.findIndex((s) => s.id === currentId);
  if (at < 0) return [];
  const ids = new Set<number>();
  const next = samples[at + 1];
  const prev = samples[at - 1];
  if (next) ids.add(next.id);
  if (prev) ids.add(prev.id);
  const unreviewed = nextUnreviewed(samples, currentId, 1);
  if (unreviewed !== undefined) ids.add(unreviewed);
  return samples.filter((s) => ids.has(s.id));
}

const encodeNeighbours = (samples: CorpusSample[]): string =>
  samples.map((s) => `${s.id}:${s.representative_exposure_id ?? ""}`).join(",");

const decodeNeighbours = (key: string): Neighbour[] =>
  key === ""
    ? []
    : key.split(",").map((part) => {
        const [id, rep] = part.split(":");
        return { id: Number(id), rep: rep ? Number(rep) : undefined };
      });

export interface PrefetchTarget {
  samples: readonly CorpusSample[];
  currentId: number | undefined;
  exposures: readonly Exposure[];
  shownId: number | undefined;
  settled: boolean;
}

export function usePrefetchNeighbours({ samples, currentId, exposures, shownId, settled }: PrefetchTarget): void {
  const qc = useQueryClient();
  const neighbourKey = encodeNeighbours(neighbourSamples(samples, currentId));
  const framesKey = exposures.flatMap((e) => imageUrl(e, false) ?? []).join("\n");
  const shown = exposures.find((e) => e.id === shownId);
  const shownUrl = useRef<string | undefined>(undefined);
  shownUrl.current = (shown && imageUrl(shown, false)) ?? undefined;

  useEffect(() => {
    if (currentId === undefined || !settled) return;
    const run = new PrefetchRun(qc);
    const timer = setTimeout(() => {
      const first = shownUrl.current;
      const frames = framesKey ? framesKey.split("\n") : [];
      void prefetchAround(
        run,
        first ? [first] : [],
        frames.filter((url) => url !== first),
        decodeNeighbours(neighbourKey),
      );
    }, PREFETCH_DWELL_MS);
    return () => {
      clearTimeout(timer);
      run.cancel();
    };
  }, [qc, currentId, settled, neighbourKey, framesKey]);
}
