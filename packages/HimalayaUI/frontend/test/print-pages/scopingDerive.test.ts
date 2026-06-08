import { describe, it, expect } from "vitest";
import {
  buildFootState,
  canScopeBuild,
  buildScopePayload,
  toPreviewSegments,
} from "../../src/print/pages/scopingDerive";
import type { OrderingRow } from "../../src/lib/scoping/proposeOrdering";
import type { PhaseRead } from "../../src/lib/scoping/dominantPhase";

const row = (over: Partial<OrderingRow>): OrderingRow => ({
  sampleId: 1,
  sampleName: "s",
  value: "1 : 0.5",
  flagged: false,
  include: true,
  ...over,
});

describe("buildFootState", () => {
  it("is ready when no flags, with the member count", () => {
    expect(buildFootState(0, 6)).toEqual({
      kind: "ready",
      text: "All 6 values confirmed — ready to build",
    });
  });
  it("warns with singular/plural flag wording", () => {
    expect(buildFootState(1, 6)).toEqual({
      kind: "warn",
      text: "1 value to check before you can build",
    });
    expect(buildFootState(2, 6).text).toBe("2 values to check before you can build");
  });
});

describe("canScopeBuild", () => {
  it("is false without an ordering key", () => {
    expect(canScopeBuild([row({})], undefined)).toBe(false);
  });
  it("is false with no included members", () => {
    expect(canScopeBuild([row({ include: false })], "ratio")).toBe(false);
  });
  it("is false when an included member is flagged", () => {
    expect(canScopeBuild([row({ flagged: true })], "ratio")).toBe(false);
  });
  it("is true when all included members are unflagged (flagged loose ignored)", () => {
    expect(
      canScopeBuild([row({ sampleId: 1 }), row({ sampleId: 2, include: false, flagged: true })], "ratio"),
    ).toBe(true);
  });
});

describe("buildScopePayload", () => {
  it("writes only included, unflagged members as {sampleId, value}", () => {
    const rows = [
      row({ sampleId: 1, value: "1 : 0" }),
      row({ sampleId: 2, flagged: true }),
      row({ sampleId: 3, include: false }),
    ];
    expect(buildScopePayload(rows)).toEqual([{ sampleId: 1, value: "1 : 0" }]);
  });
});

describe("toPreviewSegments", () => {
  it("omits coexistWith when there is no second phase", () => {
    const reads: PhaseRead[] = [{ dominant: "Pn3m", coexist: null }];
    expect(toPreviewSegments(reads)).toEqual([{ phase: "Pn3m" }]);
  });
  it("wraps the coexist partner in an array", () => {
    const reads: PhaseRead[] = [{ dominant: "Pn3m", coexist: "Lamellar" }];
    expect(toPreviewSegments(reads)).toEqual([{ phase: "Pn3m", coexistWith: ["Lamellar"] }]);
  });
  it("passes a null dominant through as an unindexed cell", () => {
    expect(toPreviewSegments([{ dominant: null, coexist: null }])).toEqual([{ phase: null }]);
  });
});
