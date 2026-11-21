import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
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
import { neighbourSamples, prefetchSample } from "../../src/workspace/prefetch";

function sample(id: number, over: Partial<CorpusSample> = {}): CorpusSample {
  return { id, experiment_id: 1, name: `S${id}`, notes: null, tags: [], q_units: "A^-1", ...over };
}

function exposure(id: number, over: Partial<Exposure> = {}): Exposure {
  return {
    id, sample_id: 1, filename: `f${id}.tif`, kind: "file", selected: false, status: null,
    image_path: "x.tif", image_version: "v1", tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

const requested: string[] = [];

class FakeImage {
  fetchPriority = "auto";
  decoding = "auto";
  onload: (() => void) | null = null;
  onerror: (() => void) | null = null;
  set src(url: string) {
    requested.push(url);
  }
}

beforeEach(() => {
  requested.length = 0;
  vi.stubGlobal("Image", FakeImage);
  vi.mocked(api.listExposures).mockResolvedValue([exposure(11), exposure(12, { status: "rejected" })]);
  vi.mocked(api.getTrace).mockResolvedValue({ q: [1, 2], I: [3, 4], sigma: [0, 0], frames: [] } as never);
  vi.mocked(api.listPeaks).mockResolvedValue([]);
  vi.mocked(api.listIndices).mockResolvedValue([]);
  vi.mocked(api.getAssignment).mockResolvedValue({ exposure_id: 12, state: "null", members: [] });
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.clearAllMocks();
});

describe("neighbourSamples", () => {
  const list = [
    sample(1, { assignment_state: "indexed", phase: "Pn3m" }),
    sample(2),
    sample(3, { assignment_state: "indexed", phase: "Im3m" }),
    sample(4, { assignment_state: "form_factor" }),
    sample(5),
  ];

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

describe("prefetchSample", () => {
  it("fills the caches the workspace reads, without waiting for the exposure list to pick the frame", async () => {
    const qc = new QueryClient();
    let releaseExposures: () => void = () => {};
    vi.mocked(api.listExposures).mockReturnValue(
      new Promise((resolve) => {
        releaseExposures = () => resolve([exposure(11)]);
      }),
    );
    const done = prefetchSample(qc, sample(1, { representative_exposure_id: 11 }), { fullImages: false });
    await Promise.resolve();
    expect(api.getTrace).toHaveBeenCalledWith(11);
    expect(api.listPeaks).toHaveBeenCalledWith(11);
    expect(api.listIndices).toHaveBeenCalledWith(11);
    expect(api.getAssignment).toHaveBeenCalledWith(11);
    releaseExposures();
    await done;
    expect(qc.getQueryData(queryKeys.exposures(1))).toEqual([exposure(11)]);
    expect(qc.getQueryData(queryKeys.trace(11))).toBeDefined();
    expect(qc.getQueryData(queryKeys.assignment(11))).toBeDefined();
  });

  it("preloads every frame's thumbnail, and full images only when asked", async () => {
    const qc = new QueryClient();
    await prefetchSample(qc, sample(1, { representative_exposure_id: 11 }), { fullImages: false });
    expect(requested).toEqual(["/api/exposures/11/image?thumb=1&v=v1", "/api/exposures/12/image?thumb=1&v=v1"]);

    requested.length = 0;
    await prefetchSample(new QueryClient(), sample(1, { representative_exposure_id: 11 }), { fullImages: true });
    expect(requested).toContain("/api/exposures/11/image?v=v1");
    expect(requested).toContain("/api/exposures/12/image?v=v1");
  });

  it("does not refetch data that is already cached and fresh", async () => {
    const qc = new QueryClient({ defaultOptions: { queries: { staleTime: 30_000 } } });
    await prefetchSample(qc, sample(1, { representative_exposure_id: 11 }), { fullImages: false });
    await prefetchSample(qc, sample(1, { representative_exposure_id: 11 }), { fullImages: false });
    expect(api.listExposures).toHaveBeenCalledTimes(1);
    expect(api.getTrace).toHaveBeenCalledTimes(1);
  });

  it("still caches the exposure list for a sample with no kept frame", async () => {
    const qc = new QueryClient();
    await prefetchSample(qc, sample(1, { representative_exposure_id: null }), { fullImages: false });
    expect(api.getTrace).not.toHaveBeenCalled();
    expect(qc.getQueryData(queryKeys.exposures(1))).toHaveLength(2);
  });
});
