import { describe, it, expect } from "vitest";
import {
  sortSampleRows,
  nextSortState,
  statusRank,
  isSortableKey,
  type SortAccessors,
  type SortState,
} from "../../../src/lib/sample/sortSamples";

interface Row {
  name: string;
  exposures: number;
  kept: number;
  phase: string | null;
}

const acc: SortAccessors<Row> = {
  name: (r) => r.name,
  exposures: (r) => r.exposures,
  kept: (r) => r.kept,
  phase: (r) => r.phase,
};

// Ingest order matters: "Sample 10" precedes "Sample 2" so a naive string sort
// (which would put 10 < 2) is distinguishable from a numeric-aware one.
const rows: Row[] = [
  { name: "Sample 10", exposures: 2, kept: 1, phase: "Pn3m" }, // 0
  { name: "Sample 2", exposures: 5, kept: 5, phase: null }, // 1
  { name: "Sample 1", exposures: 5, kept: 0, phase: "Im3m" }, // 2
  { name: "Sample 3", exposures: 1, kept: 1, phase: null }, // 3
];

const names = (rs: Row[]) => rs.map((r) => r.name);

describe("sortSampleRows", () => {
  it("a null key returns ingest order (a copy, input untouched)", () => {
    const out = sortSampleRows(rows, { key: null, dir: "asc" }, acc);
    expect(names(out)).toEqual(["Sample 10", "Sample 2", "Sample 1", "Sample 3"]);
    expect(out).not.toBe(rows);
  });

  it("Sample ascending is numeric-aware: 'Sample 2' < 'Sample 10'", () => {
    const out = sortSampleRows(rows, { key: "sample", dir: "asc" }, acc);
    expect(names(out)).toEqual(["Sample 1", "Sample 2", "Sample 3", "Sample 10"]);
  });

  it("Sample descending reverses the numeric-aware order", () => {
    const out = sortSampleRows(rows, { key: "sample", dir: "desc" }, acc);
    expect(names(out)).toEqual(["Sample 10", "Sample 3", "Sample 2", "Sample 1"]);
  });

  it("Exposures ascending sorts by frame count, ingest order breaks ties", () => {
    const out = sortSampleRows(rows, { key: "exposures", dir: "asc" }, acc);
    // counts: 2,5,5,1 → 1(Sample3), 2(Sample10), then the two 5s in ingest order
    // (Sample 2 before Sample 1).
    expect(names(out)).toEqual(["Sample 3", "Sample 10", "Sample 2", "Sample 1"]);
  });

  it("Kept ascending sorts by kept count", () => {
    const out = sortSampleRows(rows, { key: "kept", dir: "asc" }, acc);
    // kept: 1,5,0,1 → 0(Sample1), then the two 1s in ingest order (Sample10, Sample3), then 5(Sample2)
    expect(names(out)).toEqual(["Sample 1", "Sample 10", "Sample 3", "Sample 2"]);
  });

  it("Status ascending floats not-indexed (actionable) rows to the top", () => {
    const out = sortSampleRows(rows, { key: "status", dir: "asc" }, acc);
    // not-indexed (phase null): Sample 2, Sample 3 — kept in ingest order;
    // indexed: Sample 10, Sample 1 — kept in ingest order.
    expect(names(out)).toEqual(["Sample 2", "Sample 3", "Sample 10", "Sample 1"]);
  });

  it("the tiebreaker is ingest order in BOTH directions (stable)", () => {
    // Two rows tie on exposures (5). Ascending keeps Sample 2 before Sample 1;
    // descending must ALSO keep Sample 2 before Sample 1 (tiebreaker never flips).
    const asc = sortSampleRows(rows, { key: "exposures", dir: "asc" }, acc);
    const desc = sortSampleRows(rows, { key: "exposures", dir: "desc" }, acc);
    const ascTies = names(asc).filter((n) => n === "Sample 2" || n === "Sample 1");
    const descTies = names(desc).filter((n) => n === "Sample 2" || n === "Sample 1");
    expect(ascTies).toEqual(["Sample 2", "Sample 1"]);
    expect(descTies).toEqual(["Sample 2", "Sample 1"]);
  });
});

describe("statusRank", () => {
  it("not-indexed ranks before indexed (ascending = actionable first)", () => {
    expect(statusRank(null)).toBeLessThan(statusRank("Pn3m"));
    expect(statusRank("")).toBe(statusRank(undefined));
    expect(statusRank("Im3m")).toBe(statusRank("Pn3m")); // any indexed phase ties
  });
});

describe("isSortableKey", () => {
  it("accepts the four data columns, rejects tags / checkbox / junk", () => {
    for (const k of ["sample", "exposures", "kept", "status"]) {
      expect(isSortableKey(k)).toBe(true);
    }
    expect(isSortableKey("tags")).toBe(false);
    expect(isSortableKey("nonsense")).toBe(false);
  });
});

describe("nextSortState (toggle cycle)", () => {
  const at = (s: SortState) => s;

  it("inactive → ascending", () => {
    expect(nextSortState(at({ key: null, dir: "asc" }), "sample")).toEqual({
      key: "sample",
      dir: "asc",
    });
  });

  it("ascending → descending on the same column", () => {
    expect(nextSortState(at({ key: "sample", dir: "asc" }), "sample")).toEqual({
      key: "sample",
      dir: "desc",
    });
  });

  it("descending → cleared (back to null) on the same column", () => {
    expect(nextSortState(at({ key: "sample", dir: "desc" }), "sample")).toEqual({
      key: null,
      dir: "asc",
    });
  });

  it("activating a different column moves the sort there at ascending", () => {
    expect(nextSortState(at({ key: "sample", dir: "desc" }), "kept")).toEqual({
      key: "kept",
      dir: "asc",
    });
  });
});
