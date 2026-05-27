import { describe, it, expect } from "vitest";
import { proposeOrdering } from "../src/lib/scoping/proposeOrdering";
import type { SampleTagPair, PickerSampleRow } from "../src/api";

function row(
  id: number,
  name: string,
  tags: { key: string; value: string }[] = [],
): PickerSampleRow {
  return {
    sample: {
      id,
      experiment_id: 1,
      name,
      display_name: null,
      notes: null,
      tags: tags.map((t, i) => ({ id: i + 1, key: t.key, value: t.value, source: "manual" })),
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

describe("proposeOrdering", () => {
  it("proposes the most frequent corpus tag key as the ordering variable", () => {
    const corpusTags: SampleTagPair[] = [
      { key: "ratio", value: "1:1" },
      { key: "ratio", value: "2:1" },
      { key: "lipid", value: "DOPC" },
    ];
    const samples = [
      row(10, "A", [{ key: "ratio", value: "1:1" }]),
      row(11, "B", [{ key: "ratio", value: "2:1" }]),
    ];
    const result = proposeOrdering(corpusTags, samples);
    expect(result.orderingKey).toBe("ratio");
    expect(result.rows).toEqual([
      { sampleId: 10, sampleName: "A", value: "1:1", flagged: false, include: true },
      { sampleId: 11, sampleName: "B", value: "2:1", flagged: false, include: true },
    ]);
  });

  it("flags a sample missing the ordering-variable value", () => {
    const corpusTags: SampleTagPair[] = [{ key: "ratio", value: "1:1" }];
    const samples = [
      row(10, "A", [{ key: "ratio", value: "1:1" }]),
      row(11, "B", []),
    ];
    const result = proposeOrdering(corpusTags, samples);
    expect(result.rows.find((r) => r.sampleId === 11)).toEqual(
      { sampleId: 11, sampleName: "B", value: "", flagged: true, include: true });
  });

  it("returns an empty proposal on a cold corpus (no tags) — accepted by design (#174)", () => {
    const result = proposeOrdering([], [row(10, "A")]);
    expect(result.orderingKey).toBeUndefined();
    expect(result.rows).toEqual([
      { sampleId: 10, sampleName: "A", value: "", flagged: true, include: true },
    ]);
  });
});
