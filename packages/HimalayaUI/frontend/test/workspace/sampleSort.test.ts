import { describe, it, expect } from "vitest";
import type { CorpusSample } from "../../src/api";
import { sortSamples } from "../../src/workspace/sampleSort";

function sample(id: number, name: string, measured_at: string | null, setup_id = 1): CorpusSample {
  return { id, experiment_id: 1, name, notes: null, q_units: "A^-1", setup_id, measured_at };
}

const bySetup = (s: CorpusSample): number => s.setup_id ?? 0;
const names = (xs: readonly CorpusSample[]): string[] => xs.map((s) => s.name);

describe("sortSamples", () => {
  const rescanned = sample(1, "MO_AMP_10", "2026-03-02T10:00:00");
  const first = sample(2, "MO_AMP_2", "2026-03-02T09:00:00");
  const untimed = sample(3, "buffer", null);
  const later = sample(4, "mo_amp_20", "2026-03-02T11:00:00");

  it("orders samples by when their first image was taken, untimed ones last", () => {
    expect(names(sortSamples([untimed, rescanned, later, first], "measured", bySetup))).toEqual(["MO_AMP_2", "MO_AMP_10", "mo_amp_20", "buffer"]);
  });

  it("sorts names with numbers in numeric order and ignores case", () => {
    expect(names(sortSamples([later, rescanned, untimed, first], "name", bySetup))).toEqual(["buffer", "MO_AMP_2", "MO_AMP_10", "mo_amp_20"]);
  });

  it("keeps each setup's samples together in setup order", () => {
    const other = sample(5, "AAA", "2026-03-01T00:00:00", 2);
    const rank = (s: CorpusSample): number => (s.setup_id === 2 ? 1 : 0);
    expect(names(sortSamples([other, later, first], "name", rank))).toEqual(["MO_AMP_2", "mo_amp_20", "AAA"]);
    expect(names(sortSamples([other, later, first], "measured", rank))).toEqual(["MO_AMP_2", "mo_amp_20", "AAA"]);
  });

  it("breaks ties by measurement time, then by import order", () => {
    const twinLate = sample(6, "twin", "2026-03-02T12:00:00");
    const twinEarly = sample(7, "twin", "2026-03-02T08:00:00");
    const sameTimeA = sample(9, "x", "2026-03-02T08:00:00");
    const sameTimeB = sample(8, "y", "2026-03-02T08:00:00");
    expect(sortSamples([twinLate, twinEarly], "name", bySetup).map((s) => s.id)).toEqual([7, 6]);
    expect(sortSamples([sameTimeA, sameTimeB], "measured", bySetup).map((s) => s.id)).toEqual([8, 9]);
  });
});
