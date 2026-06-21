import { describe, it, expect } from "vitest";
import {
  resolveExperimentFilter,
  UNKNOWN_EXPERIMENT_LABEL,
} from "../../src/lib/experimentFilter";
import type { CorpusSample, Experiment } from "../../src/api";

// Minimal fixtures — only the fields the resolver reads.
const EXPERIMENTS = [{ id: 1, name: "BL-19 April" }] as Experiment[];
const SAMPLES = [
  { id: 10, experiment_id: 1 },
  { id: 11, experiment_id: 2 }, // experiment 2 has no /experiments record
] as CorpusSample[];

describe("resolveExperimentFilter (the SA-F5 shared predicate)", () => {
  it("no param → no filter, never unknown", () => {
    const f = resolveExperimentFilter(null, EXPERIMENTS, SAMPLES);
    expect(f).toEqual({ id: undefined, unknown: false, name: undefined });
  });

  it("malformed param → no filter (matches the page ignoring it)", () => {
    const f = resolveExperimentFilter("abc", EXPERIMENTS, SAMPLES);
    expect(f).toEqual({ id: undefined, unknown: false, name: undefined });
  });

  it("known experiment resolves its record name", () => {
    const f = resolveExperimentFilter("1", EXPERIMENTS, SAMPLES);
    expect(f).toEqual({ id: 1, unknown: false, name: "BL-19 April" });
  });

  it("an id vouched for only by samples is real, with the fallback name", () => {
    const f = resolveExperimentFilter("2", EXPERIMENTS, SAMPLES);
    // Vouched-but-unlisted: real, but there is no NAME — fallback labels are
    // each consumer's own register, not the resolver's.
    expect(f).toEqual({ id: 2, unknown: false, name: undefined });
  });

  it("an id neither list knows is unknown once both are loaded", () => {
    const f = resolveExperimentFilter("99", EXPERIMENTS, SAMPLES);
    expect(f).toEqual({ id: 99, unknown: true, name: undefined });
  });

  it("never calls a filter unknown while either list is still loading", () => {
    expect(resolveExperimentFilter("99", undefined, SAMPLES).unknown).toBe(false);
    expect(resolveExperimentFilter("99", EXPERIMENTS, undefined).unknown).toBe(false);
  });

  it("exports the single honest label both surfaces render", () => {
    expect(UNKNOWN_EXPERIMENT_LABEL).toBe("Unknown experiment");
  });
});
