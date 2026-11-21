import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { createElement, type ReactNode } from "react";
import { renderHook } from "@testing-library/react";
import { QueryClient, QueryClientProvider, useQuery } from "@tanstack/react-query";
import type { CorpusSample, Exposure } from "../../src/api";

vi.mock("../../src/api", () => ({
  listExposures: vi.fn(),
  getTrace: vi.fn(),
  listPeaks: vi.fn(),
  listIndices: vi.fn(),
  getAssignment: vi.fn(),
}));

import * as api from "../../src/api";
import { queryKeys } from "../../src/queries";
import { neighbourSamples, PREFETCH_DWELL_MS, usePrefetchNeighbours, type PrefetchTarget } from "../../src/workspace/prefetch";

function sample(id: number, over: Partial<CorpusSample> = {}): CorpusSample {
  return { id, experiment_id: 1, name: `S${id}`, notes: null, q_units: "A^-1", representative_exposure_id: id * 10 + 1, ...over };
}

function exposure(id: number, sampleId: number): Exposure {
  return {
    id, sample_id: sampleId, filename: `f${id}.tif`, kind: "file", selected: false, status: null, exposure_time: 1,
    image_path: "x.tif", image_version: "v1", sources: [], trace_hash: null, analysis_inputs_hash: null,
    counting_time: null, transmitted_counts: null, incident_monitor: null,
  };
}

const framesOf = (sampleId: number): Exposure[] => [1, 2, 3].map((n) => exposure(sampleId * 10 + n, sampleId));
const full = (id: number): string => `/api/exposures/${id}/image?v=v1`;
const thumb = (id: number): string => `/api/exposures/${id}/image?thumb=1&v=v1`;

const images: FakeImage[] = [];
const requested = (): string[] => images.map((i) => i.requested);

class FakeImage {
  decoding = "auto";
  onload: (() => void) | null = null;
  onerror: (() => void) | null = null;
  requested = "";
  current = "";
  constructor() {
    images.push(this);
  }
  set src(url: string) {
    if (url) this.requested = url;
    this.current = url;
  }
  get src(): string {
    return this.current;
  }
}

async function loadAll(): Promise<void> {
  for (let round = 0; round < 20; round++) {
    const pending = images.filter((i) => i.current && i.onload);
    if (pending.length === 0) return;
    for (const img of pending) {
      const done = img.onload!;
      img.onload = null;
      done();
    }
    await vi.advanceTimersByTimeAsync(0);
  }
}

const list = [
  sample(1, { assignment_state: "indexed", phase: "Pn3m" }),
  sample(2),
  sample(3),
  sample(4, { assignment_state: "indexed", phase: "Im3m" }),
  sample(5),
];

function target(currentId: number, over: Partial<PrefetchTarget> = {}): PrefetchTarget {
  return { samples: list, currentId, exposures: framesOf(currentId), shownId: currentId * 10 + 1, settled: true, ...over };
}

function setup(initial: PrefetchTarget, qc = new QueryClient({ defaultOptions: { queries: { retry: false } } })) {
  const wrapper = ({ children }: { children: ReactNode }) => createElement(QueryClientProvider, { client: qc }, children);
  const hook = renderHook((t: PrefetchTarget) => usePrefetchNeighbours(t), { initialProps: initial, wrapper });
  return { qc, wrapper, ...hook };
}

beforeEach(() => {
  vi.useFakeTimers();
  images.length = 0;
  vi.stubGlobal("Image", FakeImage);
  vi.mocked(api.listExposures).mockImplementation(async (id) => framesOf(id));
  vi.mocked(api.getTrace).mockResolvedValue({ q: [1, 2], I: [3, 4], sigma: [0, 0] } as never);
  vi.mocked(api.listPeaks).mockResolvedValue([]);
  vi.mocked(api.listIndices).mockResolvedValue([]);
  vi.mocked(api.getAssignment).mockImplementation(async (id) => ({ exposure_id: id, state: "null", members: [] }));
});

afterEach(() => {
  vi.useRealTimers();
  vi.unstubAllGlobals();
  vi.clearAllMocks();
});

describe("neighbourSamples", () => {
  it("returns the previous, the next and the next unreviewed sample", () => {
    expect(neighbourSamples(list, 3).map((s) => s.id)).toEqual([2, 4, 5]);
  });

  it("does not repeat a sample that is both next and next unreviewed", () => {
    expect(neighbourSamples(list, 1).map((s) => s.id)).toEqual([2]);
  });

  it("wraps the next unreviewed sample around the end of the list", () => {
    expect(neighbourSamples(list, 5).map((s) => s.id)).toEqual([2, 4]);
  });

  it("is empty when the current sample is not in the list", () => {
    expect(neighbourSamples(list, 99)).toEqual([]);
  });
});

describe("usePrefetchNeighbours", () => {
  it("does nothing until the open sample's own data has settled", async () => {
    const { rerender } = setup(target(3, { settled: false }));
    await vi.advanceTimersByTimeAsync(1000);
    expect(api.listExposures).not.toHaveBeenCalled();
    expect(api.getTrace).not.toHaveBeenCalled();
    expect(images).toHaveLength(0);

    rerender(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    expect(api.listExposures).toHaveBeenCalled();
  });

  it("waits for the user to stay on a sample before fetching anything", async () => {
    setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS - 10);
    expect(api.listExposures).not.toHaveBeenCalled();
    expect(images).toHaveLength(0);
    await vi.advanceTimersByTimeAsync(10);
    expect(api.listExposures).toHaveBeenCalled();
  });

  it("loads the shown frame's full image and the neighbours' data first, then the other frames one at a time", async () => {
    const { qc } = setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);

    expect(requested()[0]).toBe(full(31));
    expect(requested()).not.toContain(full(32));
    expect(vi.mocked(api.listExposures).mock.calls.map((c) => c[0]).sort()).toEqual([2, 4, 5]);
    expect(vi.mocked(api.getTrace).mock.calls.map((c) => c[0]).sort()).toEqual([21, 41, 51]);
    expect(qc.getQueryData(queryKeys.trace(41))).toBeDefined();
    expect(qc.getQueryData(queryKeys.assignment(51))).toBeDefined();
    expect(requested()).toEqual(expect.arrayContaining([thumb(21), thumb(22), thumb(43), thumb(53)]));

    const before32 = requested().length;
    await loadAll();
    const after = requested().slice(before32);
    expect(after.slice(0, 2)).toEqual([full(32), full(33)]);
    expect(after.slice(2).sort()).toEqual([full(21), full(41), full(51)]);
    expect(requested()).not.toContain(thumb(31));
    expect(requested()).not.toContain(full(22));
  });

  it("starts the next frame only after the previous one has loaded", async () => {
    setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    for (const img of images.filter((i) => i.requested !== full(31))) img.onload?.();
    await vi.advanceTimersByTimeAsync(0);
    expect(requested()).not.toContain(full(32));
    images.find((i) => i.requested === full(31))!.onload!();
    await vi.advanceTimersByTimeAsync(0);
    expect(requested()).toContain(full(32));
    expect(requested()).not.toContain(full(33));
  });

  it("skips a sample the user passes through quickly", async () => {
    const { rerender } = setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS - 50);
    rerender(target(4));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    expect(api.listExposures).not.toHaveBeenCalledWith(2, expect.anything());
    expect(requested()).not.toContain(full(31));
    expect(requested()[0]).toBe(full(41));
  });

  it("cancels in-flight images and unobserved queries on a move, but not a query the next sample now reads", async () => {
    const signals = new Map<number, AbortSignal>();
    vi.mocked(api.getTrace).mockImplementation((id, opts) => {
      signals.set(id, opts!.signal!);
      return new Promise(() => {});
    });
    const { qc, wrapper, rerender } = setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    const shown = images.find((i) => i.requested === full(31))!;
    expect(shown.current).toBe(full(31));

    renderHook(() => useQuery({ queryKey: queryKeys.trace(41), queryFn: () => new Promise(() => {}) }), { wrapper });
    rerender(target(4));
    await vi.advanceTimersByTimeAsync(0);

    expect(shown.current).toBe("");
    expect(signals.get(21)!.aborted).toBe(true);
    expect(signals.get(51)!.aborted).toBe(true);
    expect(signals.get(41)!.aborted).toBe(false);
    expect(qc.getQueryState(queryKeys.trace(41))?.fetchStatus).toBe("fetching");
  });

  it("keeps loading when the user steps to another frame of the same sample", async () => {
    const { rerender } = setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    const shown = images.find((i) => i.requested === full(31))!;
    const calls = vi.mocked(api.getTrace).mock.calls.length;
    rerender(target(3, { shownId: 32 }));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    expect(shown.current).toBe(full(31));
    expect(api.getTrace).toHaveBeenCalledTimes(calls);
  });

  it("does not restart when the sample list is refetched with the same neighbours", async () => {
    const { rerender } = setup(target(3));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    const count = images.length;
    rerender(target(3, { samples: list.map((s) => ({ ...s })), exposures: framesOf(3) }));
    await vi.advanceTimersByTimeAsync(PREFETCH_DWELL_MS);
    expect(images).toHaveLength(count);
    expect(api.listExposures).toHaveBeenCalledTimes(3);
  });
});
