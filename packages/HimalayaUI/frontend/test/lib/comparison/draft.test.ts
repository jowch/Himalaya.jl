import { describe, it, expect } from "vitest";
import { emptyDraft } from "../../../src/lib/comparison/draft";

// I5.3 (#184): the `fromComparison` projection cases were removed with the
// Compare-only draftFactories (orphaned after the Compare page retired). The
// kept `emptyDraft` seed-state test remains — it covers the KEPT draft.ts
// module, still consumed by the series-builder render core.
describe("ActiveDraft includes view choices — Compare UX C-4", () => {
  it("emptyDraft seeds view choices as undefined", () => {
    const d = emptyDraft();
    expect(d.viewGroupingMode).toBeUndefined();
    expect(d.viewShowPeakTicks).toBeUndefined();
    expect(d.viewShowPeakLabels).toBeUndefined();
  });
});
