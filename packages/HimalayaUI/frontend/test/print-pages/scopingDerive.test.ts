import { describe, it, expect } from "vitest";
import {
  buildFootState,
  canScopeBuild,
  buildScopePayload,
  toPreviewSegments,
  isKept,
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
  it("is ready with the kept count when nothing is skipped", () => {
    expect(buildFootState(6, 0)).toEqual({
      kind: "ready",
      text: "6 values ready to commit",
    });
  });
  it("annotates skipped members on the ready line", () => {
    expect(buildFootState(5, 1)).toEqual({
      kind: "ready",
      text: "5 values ready to commit · 1 skipped",
    });
  });
  it("uses singular wording for a single kept value", () => {
    expect(buildFootState(1, 0).text).toBe("1 value ready to commit");
    expect(buildFootState(1, 2).text).toBe("1 value ready to commit · 2 skipped");
  });
  it("warns when nothing is kept (every member skipped)", () => {
    expect(buildFootState(0, 3)).toEqual({
      kind: "warn",
      text: "Keep at least one value to build",
    });
  });
});

describe("isKept", () => {
  it("keeps an included, non-skipped row with a value", () => {
    expect(isKept(row({}))).toBe(true);
  });
  it("excludes a skipped (flagged) row", () => {
    expect(isKept(row({ flagged: true }))).toBe(false);
  });
  it("excludes a non-included row", () => {
    expect(isKept(row({ include: false }))).toBe(false);
  });
  it("excludes an empty-value row (loose matches never reach the write)", () => {
    expect(isKept(row({ value: "" }))).toBe(false);
  });
  it("IS the membership predicate behind buildScopePayload (single definition)", () => {
    const rows = [
      row({ sampleId: 1 }),
      row({ sampleId: 2, flagged: true }),
      row({ sampleId: 3, include: false }),
      row({ sampleId: 4, value: "" }),
    ];
    expect(buildScopePayload(rows).map((t) => t.sampleId)).toEqual(
      rows.filter(isKept).map((r) => r.sampleId),
    );
  });
});

describe("canScopeBuild", () => {
  it("is false without an ordering key", () => {
    expect(canScopeBuild([row({})], undefined)).toBe(false);
  });
  it("is false with no included members", () => {
    expect(canScopeBuild([row({ include: false })], "ratio")).toBe(false);
  });
  it("is false when an empty-value row is the only included member", () => {
    expect(canScopeBuild([row({ value: "" })], "ratio")).toBe(false);
  });
  it("is true when another member is kept even though one is skipped", () => {
    expect(
      canScopeBuild([row({ sampleId: 1 }), row({ sampleId: 2, flagged: true })], "ratio"),
    ).toBe(true);
  });
  it("is false when every member is skipped (flagged), never blocking elsewhere", () => {
    expect(
      canScopeBuild(
        [row({ sampleId: 1, flagged: true }), row({ sampleId: 2, flagged: true })],
        "ratio",
      ),
    ).toBe(false);
  });
});

describe("buildScopePayload", () => {
  it("writes only included, non-skipped members WITH a value as {sampleId, value}", () => {
    const rows = [
      row({ sampleId: 1, value: "1 : 0" }),
      row({ sampleId: 2, flagged: true }),
      row({ sampleId: 3, include: false }),
      row({ sampleId: 4, value: "" }),
    ];
    expect(buildScopePayload(rows)).toEqual([{ sampleId: 1, value: "1 : 0" }]);
  });
  it("never emits an empty-value entry even if flagged is false", () => {
    expect(buildScopePayload([row({ sampleId: 9, value: "", flagged: false })])).toEqual([]);
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
