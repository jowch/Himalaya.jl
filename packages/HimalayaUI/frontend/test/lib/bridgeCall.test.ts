import { describe, expect, it } from "vitest";
import { callsConflict, keepChoices } from "../../src/lib/bridgeCall";

const s = (id: number, state: "indexed" | "form_factor" | "null", phases: string[] = []) =>
  ({ id, name: `S${id}`, assignment_state: state, call_phases: phases });

describe("bridge calls", () => {
  it("two calls conflict only when both samples are called and the calls differ", () => {
    expect(callsConflict(s(1, "indexed", ["Pn3m"]), s(2, "indexed", ["Im3m"]))).toBe(true);
    expect(callsConflict(s(1, "indexed", ["Pn3m", "Im3m"]), s(2, "indexed", ["Im3m", "Pn3m"]))).toBe(false);
    expect(callsConflict(s(1, "indexed", ["Pn3m"]), s(2, "indexed"))).toBe(false);
    expect(callsConflict(s(1, "indexed", ["Pn3m"]), s(2, "null"))).toBe(false);
    expect(callsConflict(s(1, "form_factor"), s(2, "indexed", ["Pn3m"]))).toBe(true);
    expect(callsConflict(s(1, "form_factor"), s(2, "form_factor"))).toBe(false);
  });

  it("offers each sample's call, and both only when both are phase calls", () => {
    expect(keepChoices(s(1, "indexed", ["Pn3m"]), s(2, "null"))).toBeUndefined();
    expect(keepChoices(s(1, "indexed", ["Pn3m"]), s(2, "indexed", ["Im3m"]))).toEqual([
      { value: "1", label: "S1: Pn3m", keep: 1 },
      { value: "2", label: "S2: Im3m", keep: 2 },
      { value: "both", label: "Both: Im3m + Pn3m", keep: "both" },
    ]);
    expect(keepChoices(s(1, "form_factor"), s(2, "indexed", ["Im3m"]))?.map((c) => c.label)).toEqual(["S1: form factor", "S2: Im3m"]);
  });
});
