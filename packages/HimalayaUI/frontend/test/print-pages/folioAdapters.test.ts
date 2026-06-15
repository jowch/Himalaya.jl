import { describe, it, expect } from "vitest";
import { membersToSegments, toCardChrome, formatEdited, stableFigNumbers } from "../../src/print/pages/folioAdapters";
import type { SeriesMember, SeriesSummary } from "../../src/api";

function member(id: number, snap: SeriesMember["snapshot"]): SeriesMember {
  return { id, series_id: 1, exposure_id: id, display_order: id, band_height: 1, y_offset: 0,
    normalization: "none", color_override: null, label_override: null, q_window_min: null,
    q_window_max: null, peak_display: null, snapshot: snap, is_stale: false, created_by: null, created_at: null };
}
function summary(over: Partial<SeriesSummary>): SeriesSummary {
  return { id: 1, title: "LL37 ratio series", description: null, content_hash: "abc", created_by: null,
    created_at: null, updated_at: "2026-06-04T00:00:00Z", forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-06-04T00:00:00Z", author_username: "alice", member_count: 4,
    member_phases: ["Ia3d", "Im3m"], member_phase_count: 2, has_stale_members: false,
    ordering_variable: "lipid ratio", spans_experiments: false, experiment_name: null, ...over } as SeriesSummary;
}

describe("membersToSegments", () => {
  it("derives {phase, coexistWith, state} consistent with toWaterfallRows", () => {
    const segs = membersToSegments([
      member(1, { effective_peaks: [], confirmed_index: { phase: "Ia3d" } as any, analysis_inputs_hash: "h",
        confirmed_phases: [{ phase: "Ia3d" } as any, { phase: "Im3m" } as any] }),
      member(2, { effective_peaks: [], confirmed_index: { phase: "Im3m" } as any, analysis_inputs_hash: "h" }),
      member(3, { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h", assignment_state: "form_factor" }),
      member(4, null),
    ]);
    expect(segs[0]).toEqual({ phase: "Ia3d", coexistWith: ["Im3m"] });
    expect(segs[1]).toEqual({ phase: "Im3m" });
    expect(segs[2]).toEqual({ phase: null, state: "form_factor" });
    expect(segs[3]).toEqual({ phase: null });
  });
});

describe("toCardChrome", () => {
  it("derives saved-card chrome from the summary (position figure, no new-match pill)", () => {
    const c = toCardChrome(summary({}), 2, new Date("2026-06-06T00:00:00Z"));
    expect(c.figLabel).toBe("Fig. 2");
    expect(c.draft).toBe(false);
    expect(c.title).toBe("LL37 ratio series");
    expect(c.sampleCount).toBe(4);
    expect(c.variable).toBe("lipid ratio");
    expect(c.author).toBe("alice");
    expect(c.editedLabel).toBe("2 days ago");
    expect(c.notice).toBeUndefined();
    expect(c.provenance).toBeNull();
  });
  it("marks a draft and a cross-experiment series", () => {
    const c = toCardChrome(summary({ content_hash: "", spans_experiments: true, title: "" }), 1, new Date("2026-06-06T00:00:00Z"));
    expect(c.figLabel).toBe("Recipe");
    expect(c.draft).toBe(true);
    expect(c.notice).toEqual({ tone: "draft" });
    expect(c.title).toBe("Untitled series");
    expect(c.provenance).not.toBeNull();
  });
  it("names the beamtime on a single-experiment card (FOL P2-2)", () => {
    const c = toCardChrome(summary({ experiment_name: "April 2026 beamtime" }), 1, new Date("2026-06-06T00:00:00Z"));
    expect(c.provenance).toBe("April 2026 beamtime");
  });
  it("keeps the cross-experiment note when spanning, ignoring any name", () => {
    const c = toCardChrome(summary({ spans_experiments: true, experiment_name: null }), 1, new Date("2026-06-06T00:00:00Z"));
    expect(c.provenance).toBe("↔ cross-experiment · q normalized");
  });
});

describe("stableFigNumbers (FOL-FIGNUM)", () => {
  it("numbers committed series 1..N by id (creation order)", () => {
    const map = stableFigNumbers([
      summary({ id: 7 }),
      summary({ id: 2 }),
      summary({ id: 5 }),
    ]);
    expect(map.get(2)).toBe(1);
    expect(map.get(5)).toBe(2);
    expect(map.get(7)).toBe(3);
  });

  it("drafts get no number and do not consume one", () => {
    const map = stableFigNumbers([
      summary({ id: 1 }),
      summary({ id: 2, content_hash: "" }), // draft between two committed
      summary({ id: 3 }),
    ]);
    expect(map.get(1)).toBe(1);
    expect(map.has(2)).toBe(false);
    expect(map.get(3)).toBe(2); // dense: the draft did not consume "2"
  });

  it("is independent of input (view) order", () => {
    const rows = [summary({ id: 3 }), summary({ id: 1 }), summary({ id: 9 })];
    const shuffled = [rows[1]!, rows[2]!, rows[0]!];
    expect(stableFigNumbers(rows)).toEqual(stableFigNumbers(shuffled));
  });
});

describe("formatEdited", () => {
  const now = new Date("2026-06-10T00:00:00Z");
  it("formats relative times", () => {
    expect(formatEdited(null, now)).toBe("recently");
    expect(formatEdited("2026-06-10T00:00:00Z", now)).toBe("just now");
    expect(formatEdited("2026-06-09T00:00:00Z", now)).toBe("yesterday");
    expect(formatEdited("2026-06-06T00:00:00Z", now)).toBe("4 days ago");
    expect(formatEdited("2026-06-06 00:00:00", now)).toBe("4 days ago"); // SQLite space form, no Z
    expect(formatEdited("2026-06-01T00:00:00Z", now)).toBe("1 week ago"); // 9 days
    expect(formatEdited("2026-05-27T00:00:00Z", now)).toBe("2 weeks ago");
    expect(formatEdited("2026-05-10T00:00:00Z", now)).toBe("4 weeks ago"); // 31 days
  });
});
