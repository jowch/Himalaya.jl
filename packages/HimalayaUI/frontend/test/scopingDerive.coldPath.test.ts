// test/scopingDerive.coldPath.test.ts
import { describe, it, expect } from "vitest";
import {
  buildColdAssignRows, canColdBuild, buildColdScopePayload,
} from "../src/print/pages/scopingDerive";

const SEED = [
  { sampleId: 10, sampleName: "A" },
  { sampleId: 11, sampleName: "B" },
];

describe("buildColdAssignRows", () => {
  it("creates one row per seeded sample, all values empty", () => {
    const rows = buildColdAssignRows(SEED);
    expect(rows).toHaveLength(2);
    expect(rows[0]).toMatchObject({ sampleId: 10, sampleName: "A", value: "" });
    expect(rows[1]).toMatchObject({ sampleId: 11, sampleName: "B", value: "" });
  });
});

describe("canColdBuild", () => {
  it("returns false when key is empty", () => {
    expect(canColdBuild("", [{ sampleId: 10, sampleName: "A", value: "1.0" }])).toBe(false);
  });

  it("returns false when any sample has an empty value", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "" },
    ];
    expect(canColdBuild("ratio", rows)).toBe(false);
  });

  it("returns true when key and all values are non-empty", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "2.0" },
    ];
    expect(canColdBuild("ratio", rows)).toBe(true);
  });
});

describe("buildColdScopePayload", () => {
  it("returns the (sampleId, value) pairs for all rows", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "2.0" },
    ];
    expect(buildColdScopePayload(rows)).toEqual([
      { sampleId: 10, value: "1.0" },
      { sampleId: 11, value: "2.0" },
    ]);
  });
});
