import { describe, it, expect } from "vitest";
import { effectiveGroupingMode } from "../../../src/lib/comparison/effectiveGroupingMode";

describe("effectiveGroupingMode — Compare UX C-4 Step 0", () => {
  it("returns draft.viewGroupingMode when set", () => {
    const draft = { viewGroupingMode: "byPhase" } as any;
    const comparison = { view_grouping_mode: "distinct" } as any;
    expect(effectiveGroupingMode(draft, comparison)).toBe("byPhase");
  });
  it("falls back to comparison.view_grouping_mode when draft has undefined", () => {
    const draft = { viewGroupingMode: undefined } as any;
    const comparison = { view_grouping_mode: "distinct" } as any;
    expect(effectiveGroupingMode(draft, comparison)).toBe("distinct");
  });
  it("returns 'bySample' when both are null/undefined (default)", () => {
    expect(effectiveGroupingMode(null, undefined)).toBe("bySample");
    expect(effectiveGroupingMode(null, { view_grouping_mode: null } as any)).toBe("bySample");
  });
});
