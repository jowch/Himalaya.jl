import { useEffect, useRef } from "react";
import { useQueryClient, type QueryClient, type QueryKey } from "@tanstack/react-query";
import * as api from "../api";
import type { CorpusSample, Exposure } from "../api";
import { queryKeys, traceRetry } from "../queries";
import { imageUrl, nextUnreviewed } from "./model";

export const PREFETCH_DWELL_MS = 150;
export const SIBLING_WAVE = 2;
export const SIBLING_DWELL_MS = 600;

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
  private openGate: () => void = () => {};
  readonly gate = new Promise<void>((resolve) => (this.openGate = resolve));

  constructor(private readonly qc: QueryClient) {}

  get signal(): AbortSignal {
    return this.controller.signal;
  }

  query<T>(
    queryKey: QueryKey,
    fetcher: (signal: AbortSignal) => Promise<T>,
    retry?: (failures: number, error: unknown) => boolean,
  ): Promise<void> {
    if (this.signal.aborted) return Promise.resolve();
    this.keys.push(queryKey);
    return this.qc.prefetchQuery({
      queryKey,
      queryFn: ({ signal }) => fetcher(signal),
      staleTime: Infinity,
      ...(retry ? { retry } : {}),
    });
  }

  exposureData(id: number): Promise<unknown> {
    return Promise.all([
      this.query(queryKeys.trace(id), (signal) => api.getTrace(id, { signal }), traceRetry),
      this.query(queryKeys.peaks(id), (signal) => api.listPeaks(id, { signal })),
      this.query(queryKeys.indices(id), (signal) => api.listIndices(id, { signal })),
      this.query(queryKeys.assignment(id), (signal) => api.getAssignment(id, { signal })),
    ]);
  }

  async neighbourData({ id, rep }: Neighbour, thumbs: Promise<void>[]): Promise<void> {
    const data = rep === undefined ? Promise.resolve() : this.exposureData(rep);
    await this.query(queryKeys.exposures(id), (signal) => api.listExposures(id, { signal }));
    const urls = this.exposures(id).flatMap((e) => imageUrl(e, true) ?? []);
    thumbs.push(...urls.map((url) => preloadImage(url, this.signal)));
    await data;
  }

  async siblings(ids: readonly number[]): Promise<void> {
    for (let at = 0; at < ids.length; at += SIBLING_WAVE) {
      if (this.signal.aborted) return;
      await Promise.all(ids.slice(at, at + SIBLING_WAVE).map((id) => this.exposureData(id)));
    }
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

  stepped(): void {
    this.openGate();
  }

  dwell(ms: number): Promise<void> {
    return new Promise((resolve) => {
      const timer = setTimeout(resolve, ms);
      this.signal.addEventListener("abort", () => {
        clearTimeout(timer);
        resolve();
      });
    });
  }

  cancel(): void {
    this.openGate();
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

async function prefetchAround(run: PrefetchRun, shownFull: string[], otherFull: string[], siblings: number[], neighbours: Neighbour[]): Promise<void> {
  const thumbs: Promise<void>[] = [];
  const shown = run.images(shownFull);
  const ready = Promise.race([run.gate, run.dwell(SIBLING_DWELL_MS)]);
  await Promise.all(neighbours.map((n) => run.neighbourData(n, thumbs)));
  await ready;
  await run.siblings(siblings);
  await Promise.all([shown, ...thumbs]);
  await run.images(otherFull);
  await run.images(
    neighbours.flatMap((n) => {
      const rep = run.exposures(n.id).find((e) => e.id === n.rep);
      return (rep && imageUrl(rep, false)) ?? [];
    }),
  );
}

export function nearestFirst(ids: readonly number[], shownId: number | undefined): number[] {
  const at = Math.max(0, ids.indexOf(shownId ?? Number.NaN));
  return ids
    .map((id, i) => ({ id, d: i - at }))
    .filter((x) => x.d !== 0 || x.id !== shownId)
    .sort((a, b) => Math.abs(a.d) - Math.abs(b.d) || a.d - b.d)
    .map((x) => x.id);
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
  const exposureIds = exposures.map((e) => e.id).join(",");
  const shown = exposures.find((e) => e.id === shownId);
  const shownIdRef = useRef<number | undefined>(undefined);
  shownIdRef.current = shownId;
  const shownUrl = useRef<string | undefined>(undefined);
  shownUrl.current = (shown && imageUrl(shown, false)) ?? undefined;
  const runRef = useRef<{ run: PrefetchRun; shown: number | undefined } | undefined>(undefined);
  const settledFor = useRef<number | undefined>(undefined);
  if (settled && currentId !== undefined) settledFor.current = currentId;
  const ready = currentId !== undefined && settledFor.current === currentId;

  useEffect(() => {
    if (runRef.current && runRef.current.shown !== shownId) runRef.current.run.stepped();
  }, [shownId]);

  useEffect(() => {
    if (currentId === undefined || !ready) return;
    const run = new PrefetchRun(qc);
    runRef.current = { run, shown: shownIdRef.current };
    const timer = setTimeout(() => {
      const first = shownUrl.current;
      const frames = framesKey ? framesKey.split("\n") : [];
      void prefetchAround(
        run,
        first ? [first] : [],
        frames.filter((url) => url !== first),
        nearestFirst(exposureIds ? exposureIds.split(",").map(Number) : [], shownIdRef.current),
        decodeNeighbours(neighbourKey),
      );
    }, PREFETCH_DWELL_MS);
    return () => {
      clearTimeout(timer);
      run.cancel();
    };
  }, [qc, currentId, ready, neighbourKey, framesKey, exposureIds]);
}
