// test/scopingDerive.seedFilter.test.ts
import { describe, it, expect } from "vitest";
import { filterPickerBySeed } from "../src/print/pages/scopingDerive";
import type { PickerSampleRow } from "../src/api";

function row(id: number): PickerSampleRow {
  return {
    sample: { id, experiment_id: 1, name: `s${id}`, display_name: null, notes: null, tags: [] },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

describe("filterPickerBySeed", () => {
  it("returns all rows when seed is null (full-corpus path)", () => {
    const rows = [row(1), row(2), row(3)];
    expect(filterPickerBySeed(rows, null)).toEqual(rows);
  });

  it("filters to only the seeded ids when seed is provided", () => {
    const rows = [row(1), row(2), row(3)];
    const filtered = filterPickerBySeed(rows, [1, 3]);
    expect(filtered.map((r) => r.sample.id)).toEqual([1, 3]);
  });

  it("returns empty array when no seeded ids match (honest empty path)", () => {
    expect(filterPickerBySeed([row(1), row(2)], [99])).toEqual([]);
  });
});
