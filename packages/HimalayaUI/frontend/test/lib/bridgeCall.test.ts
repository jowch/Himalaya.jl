import { describe, expect, it } from "vitest";
import { bridgeSteps, callsConflict, groupKeepChoices, keepChoices } from "../../src/lib/bridgeCall";

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

  it("asks once for a group: one choice per different call, naming every sample that has it, and all of them when every call is a phase call", () => {
    const here = s(1, "indexed", ["Pn3m"]);
    expect(groupKeepChoices(here, [s(2, "null"), s(3, "indexed", ["Pn3m"])])).toBeUndefined();
    expect(groupKeepChoices(here, [s(2, "indexed", ["Im3m"]), s(3, "indexed", ["Pn3m"]), s(4, "indexed", ["Ia3d"])])).toEqual([
      { value: "1", label: "S1, S3: Pn3m", keep: 1 },
      { value: "2", label: "S2: Im3m", keep: 2 },
      { value: "4", label: "S4: Ia3d", keep: 4 },
      { value: "both", label: "All: Ia3d + Im3m + Pn3m", keep: "both" },
    ]);
    expect(groupKeepChoices(here, [s(2, "form_factor"), s(3, "null")])?.map((c) => c.label)).toEqual(["S1: Pn3m", "S2: form factor"]);
  });

  it("bridges one sample at a time, naming which call to keep only where the group and the next sample differ", () => {
    const here = s(1, "indexed", ["Pn3m"]);
    const mates = [s(2, "indexed", ["Im3m"]), s(3, "null"), s(4, "indexed", ["Pn3m"]), s(5, "indexed", ["Im3m"])];
    expect(bridgeSteps(here, mates, 1)).toEqual([{ to: 2, keep: 1 }, { to: 3 }, { to: 4 }, { to: 5, keep: 1 }]);
    expect(bridgeSteps(here, mates, 2)).toEqual([{ to: 2, keep: 2 }, { to: 3 }, { to: 4, keep: 1 }, { to: 5 }]);
    expect(bridgeSteps(here, mates, "both")).toEqual([{ to: 2, keep: "both" }, { to: 3 }, { to: 4, keep: "both" }, { to: 5, keep: "both" }]);
    expect(bridgeSteps(s(1, "null"), [s(2, "indexed", ["Im3m"]), s(3, "indexed", ["Pn3m"])], 3)).toEqual([{ to: 2 }, { to: 3, keep: 3 }]);
  });
});
