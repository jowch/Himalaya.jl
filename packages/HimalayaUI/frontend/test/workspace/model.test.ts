import { describe, it, expect } from "vitest";
import type { CorpusSample, Exposure, IndexEntry } from "../../src/api";
import { sampleStatus, nextUnreviewed, rankCandidates, imageUrl, stepExposure } from "../../src/workspace/model";

function sample(id: number, over: Partial<CorpusSample> = {}): CorpusSample {
  return { id, experiment_id: 1, name: `S${id}`, notes: null, tags: [], q_units: "A^-1", ...over };
}

function index(id: number, phase: string, score: number | null): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score, r_squared: null, lattice_d: null, ngc: null,
    status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}

function exposure(id: number, over: Partial<Exposure> = {}): Exposure {
  return {
    id, sample_id: 1, filename: `f${id}.tif`, kind: "file", selected: false, status: null,
    image_path: "x.tif", image_version: "v1", tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

describe("sampleStatus", () => {
  it("is unreviewed without a confirmed call", () => {
    expect(sampleStatus(sample(1))).toEqual({ kind: "unreviewed" });
  });

  it("is unreviewed when indexed carries no phase", () => {
    expect(sampleStatus(sample(1, { assignment_state: "indexed", phase: null }))).toEqual({ kind: "unreviewed" });
  });

  it("reports the dominant assigned phase", () => {
    expect(sampleStatus(sample(1, { assignment_state: "indexed", phase: "Pn3m" }))).toEqual({
      kind: "phase", phase: "Pn3m",
    });
  });

  it("reports form factor", () => {
    expect(sampleStatus(sample(1, { assignment_state: "form_factor" }))).toEqual({ kind: "form_factor" });
  });
});

describe("nextUnreviewed", () => {
  const list = [
    sample(1, { assignment_state: "indexed", phase: "Pn3m" }),
    sample(2),
    sample(3, { assignment_state: "form_factor" }),
    sample(4),
  ];

  it("finds the next unreviewed sample after the current one, wrapping", () => {
    expect(nextUnreviewed(list, 2, 1)).toBe(4);
    expect(nextUnreviewed(list, 4, 1)).toBe(2);
  });

  it("searches backwards", () => {
    expect(nextUnreviewed(list, 4, -1)).toBe(2);
    expect(nextUnreviewed(list, 2, -1)).toBe(4);
  });

  it("returns undefined when every sample is reviewed", () => {
    expect(nextUnreviewed([list[0]!, list[2]!], 1, 1)).toBeUndefined();
  });
});

describe("rankCandidates", () => {
  it("sorts by score, drops unscored candidates, and keeps at most nine", () => {
    const ranked = rankCandidates([
      index(1, "Im3m", 0.4), index(2, "Pn3m", null), index(3, "Ia3d", 0.9),
      ...Array.from({ length: 10 }, (_, i) => index(10 + i, "Lamellar", 0.1)),
    ]);
    expect(ranked.map((c) => c.id).slice(0, 2)).toEqual([3, 1]);
    expect(ranked).toHaveLength(9);
    expect(ranked.some((c) => c.id === 2)).toBe(false);
  });
});

describe("stepExposure", () => {
  const list = [exposure(1), exposure(2), exposure(3)];

  it("moves within the list without wrapping", () => {
    expect(stepExposure(list, 1, 1)).toBe(2);
    expect(stepExposure(list, 3, 1)).toBe(3);
    expect(stepExposure(list, 1, -1)).toBe(1);
  });
});

describe("imageUrl", () => {
  it("builds a thumbnail url with the cache-busting version", () => {
    expect(imageUrl(exposure(7), true)).toBe("/api/exposures/7/image?thumb=1&v=v1");
  });

  it("is null without an image", () => {
    expect(imageUrl(exposure(7, { image_path: null }), true)).toBeNull();
  });
});
