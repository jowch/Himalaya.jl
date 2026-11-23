import { describe, it, expect } from "vitest";
import type { CorpusSample, Exposure, IndexEntry } from "../../src/api";
import { i0BelowZero, traceIsDark, sampleStatus, nextUnreviewed, firstUnreviewed, rankCandidates, imageUrl, stepExposure, nearestKeptExposure, peakClick } from "../../src/workspace/model";
import { peakList } from "../fixtures/peaks";
import { makeExposure } from "../helpers/exposure";

function sample(id: number, over: Partial<CorpusSample> = {}): CorpusSample {
  return { id, experiment_id: 1, name: `S${id}`, notes: null, q_units: "A^-1", ...over };
}

function index(id: number, phase: string, score: number | null): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score, r_squared: null, lattice_d: null, 
    status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
  };
}

function exposure(id: number, over: Partial<Exposure> = {}): Exposure {
  return makeExposure({ id, filename: `f${id}.tif`, image_path: "x.tif", image_version: "v1", ...over });
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
      kind: "phase", phases: ["Pn3m"],
    });
  });

  it("reports every phase of a coexistence call, dominant first", () => {
    expect(sampleStatus(sample(1, { assignment_state: "indexed", phase: "Im3m", call_phases: ["Pn3m", "Im3m"] }))).toEqual({
      kind: "phase", phases: ["Im3m", "Pn3m"],
    });
  });

  it("reports form factor", () => {
    expect(sampleStatus(sample(1, { assignment_state: "form_factor" }))).toEqual({ kind: "form_factor" });
  });
});

describe("firstUnreviewed", () => {
  it("opens on the first unreviewed sample", () => {
    expect(firstUnreviewed([sample(1, { assignment_state: "indexed", phase: "Pn3m" }), sample(2), sample(3)])?.id).toBe(2);
  });

  it("falls back to the first sample when everything is reviewed", () => {
    expect(firstUnreviewed([sample(1, { assignment_state: "form_factor" }), sample(2, { assignment_state: "form_factor" })])?.id).toBe(1);
  });

  it("is undefined for an empty experiment", () => {
    expect(firstUnreviewed([])).toBeUndefined();
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
  it("sorts by score, drops unscored candidates, and keeps every scored one", () => {
    const ranked = rankCandidates([
      index(1, "Im3m", 0.4), index(2, "Pn3m", null), index(3, "Ia3d", 0.9),
      ...Array.from({ length: 10 }, (_, i) => index(10 + i, "Lamellar", 0.1)),
    ]);
    expect(ranked.map((c) => c.id).slice(0, 2)).toEqual([3, 1]);
    expect(ranked).toHaveLength(12);
    expect(ranked.some((c) => c.id === 2)).toBe(false);
  });

  it("keeps a custom index, which has no score, and lists it first", () => {
    const custom = { ...index(5, "Hexagonal", null), kind: "speculative" as const };
    const ranked = rankCandidates([index(1, "Im3m", 0.4), custom, index(3, "Ia3d", 0.9)]);
    expect(ranked.map((c) => c.id)).toEqual([5, 3, 1]);
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

describe("nearestKeptExposure", () => {
  const kept = (id: number) => exposure(id);
  const dropped = (id: number) => exposure(id, { status: "rejected" });

  it("is the next kept frame after the one just dropped", () => {
    expect(nearestKeptExposure([kept(1), kept(2), kept(3)], 2)).toBe(3);
  });

  it("falls back to the previous kept frame when the dropped one is last", () => {
    expect(nearestKeptExposure([kept(1), kept(2), kept(3)], 3)).toBe(2);
  });

  it("skips frames that are already dropped, picking the closer side", () => {
    expect(nearestKeptExposure([kept(1), dropped(2), kept(3), dropped(4), dropped(5), kept(6)], 4)).toBe(3);
    expect(nearestKeptExposure([kept(1), dropped(2), dropped(3), kept(4)], 2)).toBe(1);
  });

  it("is undefined when no other frame is kept", () => {
    expect(nearestKeptExposure([dropped(1), kept(2)], 2)).toBeUndefined();
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

describe("peakClick", () => {
  const peaks = peakList(7, { auto: [0.1, { q: 0.2, excluded: true }], manual: [0.15, 0.25] });
  it("disables an active auto peak and restores a disabled one", () => {
    expect(peakClick(peaks, { id: 1, source: "auto" })).toEqual({ kind: "exclude", excluded: true });
    expect(peakClick(peaks, { id: 2, source: "auto" })).toEqual({ kind: "exclude", excluded: false });
  });
  it("removes a manual peak that shares its id with an auto peak", () => {
    expect(peakClick(peaks, { id: 1, source: "manual" })).toEqual({ kind: "remove" });
    expect(peakClick(peaks, { id: 2, source: "manual" })).toEqual({ kind: "remove" });
  });
  it("ignores an unknown peak", () => {
    expect(peakClick(peaks, { id: 3, source: "auto" })).toBeNull();
  });
});

describe("traceIsDark", () => {
  it("is true when the median count is under one", () => {
    expect(traceIsDark({ I: [0, 0, 0.002, 0, 40] })).toBe(true);
  });

  it("is false for a dim but real frame, an empty trace and no trace", () => {
    expect(traceIsDark({ I: [8, 12, 15, 0, 0] })).toBe(false);
    expect(traceIsDark({ I: [] })).toBe(false);
    expect(traceIsDark(undefined)).toBe(false);
  });
});

describe("i0BelowZero", () => {
  const f = (incident_monitor: number | null) => ({ incident_monitor });

  it("is true at or below zero", () => {
    for (const i0 of [-0.1, -0.005, 0]) expect(i0BelowZero(f(i0))).toBe(true);
  });

  it("is false for any positive reading and an unknown one", () => {
    for (const i0 of [0.004, 0.09, 1.5]) expect(i0BelowZero(f(i0))).toBe(false);
    expect(i0BelowZero(f(null))).toBe(false);
  });
});
