import { describe, it, expect } from "vitest";
import { emptyDraft } from "../../../src/lib/comparison/draft";
import { fromComparison } from "../../../src/lib/comparison/draftFactories";
import type { Comparison } from "../../../src/api";

// Minimal QueryClient stub — fromComparison calls memberFromSaved which
// calls computeMemberSnapshot(exposure_id, qc). With no members the qc
// is never touched, so a bare object suffices.
const stubQc = {} as Parameters<typeof fromComparison>[1];

describe("ActiveDraft includes view choices — Compare UX C-4", () => {
  it("emptyDraft seeds view choices as undefined", () => {
    const d = emptyDraft();
    expect(d.viewGroupingMode).toBeUndefined();
    expect(d.viewShowPeakTicks).toBeUndefined();
    expect(d.viewShowPeakLabels).toBeUndefined();
  });

  it("fromComparison projects existing view choices", () => {
    const c: Comparison = {
      id: 1,
      title: "x",
      description: null,
      content_hash: "h",
      created_by: null,
      created_at: null,
      updated_at: null,
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: false,
      view_show_peak_labels: true,
      last_event_at: null,
      members: [],
    };
    const d = fromComparison(c, stubQc);
    expect(d.viewGroupingMode).toBe("byPhase");
    expect(d.viewShowPeakTicks).toBe(false);
    expect(d.viewShowPeakLabels).toBe(true);
  });

  it("fromComparison coerces null view choices to undefined", () => {
    const c: Comparison = {
      id: 2,
      title: "y",
      description: null,
      content_hash: "h2",
      created_by: null,
      created_at: null,
      updated_at: null,
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      view_grouping_mode: null,
      view_show_peak_ticks: null,
      view_show_peak_labels: null,
      last_event_at: null,
      members: [],
    };
    const d = fromComparison(c, stubQc);
    expect(d.viewGroupingMode).toBeUndefined();
    expect(d.viewShowPeakTicks).toBeUndefined();
    expect(d.viewShowPeakLabels).toBeUndefined();
  });
});
