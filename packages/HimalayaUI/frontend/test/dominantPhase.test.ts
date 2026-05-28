import { describe, it, expect } from "vitest";
import { dominantPhase } from "../src/lib/scoping/dominantPhase";
import type { IndexEntry } from "../src/api";

const idx = (phase: string, score: number): IndexEntry => ({
  id: 1, exposure_id: 1, phase, basis: 1, score, r_squared: null, lattice_d: null,
  ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [],
});

describe("dominantPhase", () => {
  it("picks the highest-scored phase as dominant", () => {
    expect(dominantPhase([idx("Lamellar", 0.4), idx("Pn3m", 0.9)]).dominant).toBe("Pn3m");
  });
  it("reports the second distinct phase as coexist", () => {
    const r = dominantPhase([idx("Pn3m", 0.9), idx("Lamellar", 0.6)]);
    expect(r.dominant).toBe("Pn3m");
    expect(r.coexist).toBe("Lamellar");
  });
  it("no coexist when only one phase is present", () => {
    expect(dominantPhase([idx("Pn3m", 0.9)]).coexist).toBeNull();
  });
  it("returns null dominant for an unindexed exposure", () => {
    expect(dominantPhase([]).dominant).toBeNull();
  });
});
