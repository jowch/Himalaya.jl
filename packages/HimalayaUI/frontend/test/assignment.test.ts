import { describe, it, expect } from "vitest";
import { deriveActiveIndices, assignmentState } from "../src/lib/assignment";
import type { Assignment, IndexEntry } from "../src/api";

function ix(id: number, phase = "Pn3m"): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score: 0.9, r_squared: 0.99,
    lattice_d: 100, ngc: null, status: "candidate", kind: "auto",
    inputs_hash: null, peaks: [], predicted_q: [],
  };
}

const indices = [ix(10), ix(11), ix(12)];

describe("deriveActiveIndices (Plan D-2)", () => {
  it("returns [] when assignment is undefined", () => {
    expect(deriveActiveIndices(undefined, indices)).toEqual([]);
  });

  it("returns [] for form_factor state regardless of members", () => {
    const a: Assignment = { exposure_id: 1, state: "form_factor", members: [10, 11] };
    expect(deriveActiveIndices(a, indices)).toEqual([]);
  });

  it("returns [] for null state", () => {
    const a: Assignment = { exposure_id: 1, state: "null", members: [10] };
    expect(deriveActiveIndices(a, indices)).toEqual([]);
  });

  it("filters indices by members for indexed state", () => {
    const a: Assignment = { exposure_id: 1, state: "indexed", members: [11, 12] };
    expect(deriveActiveIndices(a, indices).map((i) => i.id)).toEqual([11, 12]);
  });

  it("indexed with 0 members yields [] (call in progress)", () => {
    const a: Assignment = { exposure_id: 1, state: "indexed", members: [] };
    expect(deriveActiveIndices(a, indices)).toEqual([]);
  });

  it("ignores member ids not present in the indices list", () => {
    const a: Assignment = { exposure_id: 1, state: "indexed", members: [10, 999] };
    expect(deriveActiveIndices(a, indices).map((i) => i.id)).toEqual([10]);
  });
});

describe("assignmentState (Plan D-2)", () => {
  it("defaults to indexed when undefined", () => {
    expect(assignmentState(undefined)).toBe("indexed");
  });
  it("passes through the assignment state", () => {
    expect(assignmentState({ exposure_id: 1, state: "null", members: [] })).toBe("null");
    expect(assignmentState({ exposure_id: 1, state: "form_factor", members: [] })).toBe("form_factor");
  });
});
