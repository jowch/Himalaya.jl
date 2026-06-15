import { describe, it, expect } from "vitest";
import { deriveExperimentSiblings } from "../../../src/lib/sample/experimentSiblings";
import type { CorpusSample } from "../../../src/api";

// F5: the single shared derivation behind BOTH the topbar sample stepper and
// the `,`/`.` global shortcut on the focus route — the active sample's
// experiment-siblings in corpus order. One source, so the two can never
// disagree (alignment invariant).

function cs(id: number, experiment_id: number): CorpusSample {
  return {
    id,
    experiment_id,
    name: `smp_${id}`,
    display_name: null,
    notes: null,
    tags: [],
    q_units: "A-1",
  } as CorpusSample;
}

// Corpus order is interleaved across experiments on purpose: siblings must
// preserve corpus order while filtering to the active sample's experiment.
const CORPUS = [cs(10, 1), cs(90, 2), cs(11, 1), cs(12, 1), cs(91, 2)];

describe("deriveExperimentSiblings", () => {
  it("filters to the active sample's experiment, preserving corpus order", () => {
    const d = deriveExperimentSiblings(CORPUS, 11);
    expect(d.activeSample?.id).toBe(11);
    expect(d.siblings.map((s) => s.id)).toEqual([10, 11, 12]);
    expect(d.index).toBe(1);
  });

  it("exposes prev/next within the experiment (no wrap)", () => {
    const mid = deriveExperimentSiblings(CORPUS, 11);
    expect(mid.prev?.id).toBe(10);
    expect(mid.next?.id).toBe(12);

    const first = deriveExperimentSiblings(CORPUS, 10);
    expect(first.prev).toBeUndefined();
    expect(first.next?.id).toBe(11);

    const last = deriveExperimentSiblings(CORPUS, 12);
    expect(last.prev?.id).toBe(11);
    expect(last.next).toBeUndefined();
  });

  it("never crosses experiments: an experiment-2 sample steps among exp-2 only", () => {
    const d = deriveExperimentSiblings(CORPUS, 90);
    expect(d.siblings.map((s) => s.id)).toEqual([90, 91]);
    expect(d.prev).toBeUndefined();
    expect(d.next?.id).toBe(91);
  });

  it("cold corpus cache (undefined samples) -> empty derivation, no throw", () => {
    const d = deriveExperimentSiblings(undefined, 11);
    expect(d.activeSample).toBeUndefined();
    expect(d.siblings).toEqual([]);
    expect(d.index).toBe(-1);
    expect(d.prev).toBeUndefined();
    expect(d.next).toBeUndefined();
  });

  it("unknown active sample -> empty derivation", () => {
    const d = deriveExperimentSiblings(CORPUS, 99999);
    expect(d.activeSample).toBeUndefined();
    expect(d.siblings).toEqual([]);
    expect(d.prev).toBeUndefined();
    expect(d.next).toBeUndefined();
  });

  it("no active sample id -> empty derivation", () => {
    const d = deriveExperimentSiblings(CORPUS, undefined);
    expect(d.activeSample).toBeUndefined();
    expect(d.siblings).toEqual([]);
    expect(d.index).toBe(-1);
  });
});
