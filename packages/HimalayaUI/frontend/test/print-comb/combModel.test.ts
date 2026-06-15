import { describe, it, expect } from "vitest";
import { assembleRows, combQDomain, type CombSeries } from "../../src/print/comb/combModel";

const pn3m: CombSeries = {
  phase: "Pn3m",
  color: "var(--a)",
  teeth: [
    { q: 0.10, label: "√2", observed: true },
    { q: 0.1225, label: "√3", observed: false },
  ],
};
const im3m: CombSeries = {
  phase: "Im3m",
  color: "var(--b)",
  teeth: [{ q: 0.12, label: "√2", observed: true }],
};

describe("assembleRows", () => {
  it("orders assigned rows, then the preview row, then the leftover row", () => {
    const rows = assembleRows([pn3m, im3m], im3m, [0.2]);
    expect(rows.map((r) => r.kind)).toEqual(["assigned", "assigned", "preview", "leftover"]);
  });

  it("omits the preview row when there is no hovered series", () => {
    const rows = assembleRows([pn3m], undefined, [0.2]);
    expect(rows.some((r) => r.kind === "preview")).toBe(false);
  });

  it("omits the leftover row when there are no leftover peaks", () => {
    const rows = assembleRows([pn3m], undefined, []);
    expect(rows.some((r) => r.kind === "leftover")).toBe(false);
  });

  it("carries the leftover q-values on the leftover row", () => {
    const rows = assembleRows([pn3m], undefined, [0.2, 0.25]);
    const leftover = rows.find((r) => r.kind === "leftover");
    expect(leftover && leftover.kind === "leftover" && leftover.qs).toEqual([0.2, 0.25]);
  });
});

describe("combQDomain", () => {
  it("spans every q (teeth + leftover) with ~10% pad and ignores non-positive", () => {
    const rows = assembleRows([pn3m], undefined, [0.2, -1, 0]);
    const [lo, hi] = combQDomain(rows);
    expect(lo).toBeCloseTo(0.10 * 0.9, 6);
    expect(hi).toBeCloseTo(0.20 * 1.1, 6);
  });
});
