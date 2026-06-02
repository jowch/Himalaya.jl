import { describe, it, expect } from "vitest";
import { cursorCandidates, snapToPeakQ } from "../../src/print/waterfall/cursor";
import type { WaterfallRow } from "../../src/print/waterfall/waterfallModel";

const EMPTY = { q: [], I: [], sigma: [] };
function row(key: string, anchors: { id: number; q: number }[]): WaterfallRow {
  return {
    key, label: key, trace: EMPTY, phase: "Ia3d", state: "indexed",
    anchors: anchors.map((a) => ({ id: a.id, q: a.q, intensity: 100, phase: "Ia3d" })),
    bandHeight: 1, yOffset: 0,
  };
}
const toPx = (q: number) => q * 1000;

describe("cursor", () => {
  it("flattens anchors across all rows into {id,q} candidates", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }]), row("b", [{ id: 2, q: 0.08 }, { id: 3, q: 0.10 }])];
    expect(cursorCandidates(rows)).toEqual([
      { id: 1, q: 0.06 }, { id: 2, q: 0.08 }, { id: 3, q: 0.10 },
    ]);
  });

  it("snaps to the nearest peak q within tolerance", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }, { id: 2, q: 0.09 }])];
    expect(snapToPeakQ(62, rows, toPx, 10)).toBe(0.06);
  });

  it("returns null when no peak is within tolerance", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }])];
    expect(snapToPeakQ(100, rows, toPx, 10)).toBeNull();
  });

  it("returns null for empty rows", () => {
    expect(snapToPeakQ(60, [], toPx, 10)).toBeNull();
  });

  it("ignores optimistic peaks (id < 0)", () => {
    const rows = [row("a", [{ id: -1, q: 0.06 }])];
    expect(snapToPeakQ(60, rows, toPx, 10)).toBeNull();
  });
});
