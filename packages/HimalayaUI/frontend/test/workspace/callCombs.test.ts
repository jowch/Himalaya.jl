import { describe, it, expect } from "vitest";
import { callCombs } from "../../src/workspace/callCombs";
import { phaseColor } from "../../src/phases";
import type { IndexEntry } from "../../src/api";

const entry = (id: number, phase: string, predicted_q: number[], observed: (number | null)[] = []): IndexEntry => ({
  id,
  phase,
  predicted_q,
  basis: predicted_q[0],
  peaks: observed.flatMap((q, i) => (q === null ? [] : [{ peak_id: i + 1, ratio_position: i + 1, residual: 0, q_observed: q }])),
} as unknown as IndexEntry);
const pn3m = entry(1, "Pn3m", [0.1, 0.12]);
const im3m = entry(2, "Im3m", [0.08, 0.11]);
const hex = entry(3, "Hexagonal", [0.09]);
const candidates = [pn3m, im3m, hex];

describe("callCombs", () => {
  it("draws nothing with no call and no hover", () => {
    expect(callCombs(candidates, undefined, new Set(), null)).toEqual([]);
  });

  it("keeps the called phases drawn faded when nothing is hovered", () => {
    expect(callCombs(candidates, undefined, new Set([1, 3]), null)).toEqual([
      { id: 1, qs: [0.1, 0.12], color: phaseColor("Pn3m"), faded: true },
      { id: 3, qs: [0.09], color: phaseColor("Hexagonal"), faded: true },
    ]);
    expect(callCombs(candidates, [], new Set([1]), null)).toEqual([{ id: 1, qs: [0.1, 0.12], color: phaseColor("Pn3m"), faded: true }]);
  });

  it("draws only the hovered phases at full strength while hovering, even over a call", () => {
    expect(callCombs(candidates, [2], new Set([1]), null)).toEqual([{ id: 2, qs: [0.08, 0.11], color: phaseColor("Im3m") }]);
    expect(callCombs(candidates, [1], new Set([1]), null)).toEqual([{ id: 1, qs: [0.1, 0.12], color: phaseColor("Pn3m") }]);
  });

  it("draws only the custom draft while one is open", () => {
    expect(callCombs(candidates, [2], new Set([1]), { phase: "Ia3d", qs: [0.05] })).toEqual([{ id: -1, qs: [0.05], color: phaseColor("Ia3d") }]);
  });

  it("draws the reflections from the lattice fitted to every claimed peak, not from the first peak alone", () => {
    const ratios = [Math.sqrt(2), 2, Math.sqrt(6), Math.sqrt(8), Math.sqrt(10)].map((r) => r / Math.sqrt(2));
    const predicted = ratios.map((r) => 0.0644 * r);
    const observed = [0.0644, 0.0909, 0.1113, 0.1287, null];
    const drifting = entry(9, "Im3m", predicted, observed);
    const [comb] = callCombs([drifting], undefined, new Set([9]), null);
    const hits = observed.flatMap((q, i) => (q === null ? [] : [{ r: ratios[i]!, q }]));
    const q1 = hits.reduce((a, h) => a + h.r * h.q, 0) / hits.reduce((a, h) => a + h.r * h.r, 0);
    comb!.qs.forEach((q, i) => expect(q).toBeCloseTo(q1 * ratios[i]!, 9));
    const worst = (qs: readonly number[]): number => Math.max(...observed.map((o, i) => (o === null ? 0 : Math.abs(qs[i]! - o))));
    expect(worst(comb!.qs)).toBeLessThan(worst(predicted));
    expect(callCombs([drifting], [9], new Set(), null)[0]!.qs).toEqual(comb!.qs);
  });

  it("keeps the stored predictions when the phase claims no peak", () => {
    expect(callCombs([entry(4, "Pn3m", [0.1, 0.12])], undefined, new Set([4]), null)[0]!.qs).toEqual([0.1, 0.12]);
  });

  it("keeps a custom index's lines on the lattice the user chose, even when its claimed peaks sit off them", () => {
    const predicted = [0.1, 0.1 * Math.sqrt(1.5), 0.1 * Math.sqrt(2)];
    const custom = { ...entry(5, "Pn3m", predicted, [0.101, 0.1 * Math.sqrt(1.5) * 1.015, 0.1 * Math.sqrt(2) * 1.02]), kind: "speculative" } as IndexEntry;
    expect(callCombs([custom], undefined, new Set([5]), null)[0]!.qs).toEqual(predicted);
    expect(callCombs([custom], [5], new Set(), null)[0]!.qs).toEqual(predicted);
  });
});
