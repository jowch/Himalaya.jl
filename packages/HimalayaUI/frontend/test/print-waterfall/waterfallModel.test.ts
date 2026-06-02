import { describe, it, expect } from "vitest";
import { toWaterfallRows, waterfallQDomain } from "../../src/print/waterfall/waterfallModel";
import { realMembers, formFactorMember, unindexedMember } from "../../src/print/fixtures/realSeriesMembers";
import { realTraces } from "../../src/print/fixtures/realTraces";

describe("toWaterfallRows", () => {
  it("returns one row per member, preserving display order", () => {
    const rows = toWaterfallRows(realMembers, realTraces);
    expect(rows.length).toBe(realMembers.length);
    expect(rows.map((r) => r.key)).toEqual(realMembers.map((m) => String(m.id)));
  });

  it("derives the dominant phase from confirmed_phases / confirmed_index", () => {
    const rows = toWaterfallRows([realMembers[0]!], realTraces);
    expect(rows[0]!.phase).toBe("Ia3d"); // exp 65
  });

  it("emits an anchor per confirmed_index.peak_id, resolved to its effective-peak q", () => {
    const rows = toWaterfallRows([realMembers[0]!], realTraces);
    const ci = realMembers[0]!.snapshot!.confirmed_index!;
    expect(rows[0]!.anchors.map((a) => a.id).sort()).toEqual([...ci.peak_ids].sort());
    const byId = new Map(realMembers[0]!.snapshot!.effective_peaks.map((p) => [p.id, p.q]));
    for (const a of rows[0]!.anchors) expect(a.q).toBe(byId.get(a.id));
    const iById = new Map(realMembers[0]!.snapshot!.effective_peaks.map((p) => [p.id, p.intensity]));
    for (const a of rows[0]!.anchors) expect(a.intensity).toBe(iById.get(a.id));
  });

  it("yields zero anchors and a null phase for a form-factor member", () => {
    const rows = toWaterfallRows([formFactorMember], realTraces);
    expect(rows[0]!.anchors).toEqual([]);
    expect(rows[0]!.phase).toBeNull();
    expect(rows[0]!.state).toBe("form_factor");
  });

  it("yields zero anchors for an unindexed (null-assignment) member", () => {
    const rows = toWaterfallRows([unindexedMember], realTraces);
    expect(rows[0]!.anchors).toEqual([]);
    expect(rows[0]!.phase).toBeNull();
  });

  it("silently drops confirmed_index.peak_ids absent from effective_peaks", () => {
    const rows = toWaterfallRows([realMembers[1]!], realTraces);
    const presentIds = new Set(realMembers[1]!.snapshot!.effective_peaks.map((p) => p.id));
    for (const a of rows[0]!.anchors) expect(presentIds.has(a.id)).toBe(true);
    expect(rows[0]!.anchors.length).toBeLessThan(
      realMembers[1]!.snapshot!.confirmed_index!.peak_ids.length,
    );
  });

  it("uses an empty trace when no measured trace is found (no throw)", () => {
    const orphan = { ...realMembers[0]!, exposure_id: 999999 };
    const rows = toWaterfallRows([orphan], realTraces);
    expect(rows[0]!.trace.q).toEqual([]);
    expect(rows[0]!.trace.I).toEqual([]);
  });
});

describe("waterfallQDomain", () => {
  it("returns the padded positive q-extent across all rows", () => {
    const rows = toWaterfallRows(realMembers, realTraces);
    const [lo, hi] = waterfallQDomain(rows);
    const allQ = rows.flatMap((r) => r.trace.q).filter((q) => q > 0);
    expect(lo).toBeLessThanOrEqual(Math.min(...allQ));
    expect(hi).toBeGreaterThanOrEqual(Math.max(...allQ));
    expect(lo).toBeGreaterThan(0);
  });

  it("falls back to [0.01, 1] when there is no positive data", () => {
    expect(waterfallQDomain([])).toEqual([0.01, 1]);
  });
});
