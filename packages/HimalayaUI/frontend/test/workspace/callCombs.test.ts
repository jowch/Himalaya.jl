import { describe, it, expect } from "vitest";
import { callCombs } from "../../src/workspace/callCombs";
import { phaseColor } from "../../src/phases";
import type { IndexEntry } from "../../src/api";

const entry = (id: number, phase: string, predicted_q: number[]): IndexEntry => ({ id, phase, predicted_q } as unknown as IndexEntry);
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
});
