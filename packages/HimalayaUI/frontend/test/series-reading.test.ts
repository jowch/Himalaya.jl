/**
 * seriesReading (Plan E, Task E-4) — derive the "phases present" reading from
 * per-member snapshots: each indexed phase's variable span + lattice trend,
 * plus coexistence and form-factor lines. Generalizes to any series (no
 * hand-written X→Y narration).
 */
import { describe, it, expect } from "vitest";
import type { SeriesMember, MemberSnapshotPhase, AssignmentState } from "../src/api";
import { seriesReading } from "../src/lib/series/seriesReading";

let nextId = 1;
function member(
  phases: MemberSnapshotPhase[],
  state: AssignmentState = phases.length > 0 ? "indexed" : "form_factor",
): SeriesMember {
  const id = nextId++;
  const primary = phases[0] ?? null;
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: {
      effective_peaks: primary
        ? [{ id: id * 100, q: 0.1, intensity: 1, sharpness: 1, source: "auto" as const }]
        : [],
      confirmed_index: primary
        ? {
            id: id * 1000, phase: primary.phase, lattice_d: primary.lattice_d ?? 0,
            r_squared: 0.99, ngc: -1.5, peak_ids: [id * 100],
          }
        : null,
      confirmed_phases: phases,
      assignment_state: state,
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: null, created_at: null,
  };
}

describe("seriesReading", () => {
  it("derives each phase's variable span + lattice trend", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Pn3m", lattice_d: 205 }]),
      member([{ phase: "Pn3m", lattice_d: 202 }]),
      member([{ phase: "Pn3m", lattice_d: 198 }]),
      member([{ phase: "Pn3m", lattice_d: 195 }]),
    ];
    const vars = ["1:0", "1:0.1", "1:0.25", "1:0.5"];
    const reading = seriesReading(members, (m) => vars[m.display_order]!);
    expect(reading.phases).toHaveLength(1);
    const p = reading.phases[0]!;
    expect(p.phase).toBe("Pn3m");
    expect(p.spanLabel).toBe("1:0 → 1:0.5");
    expect(p.latticeTrend).toBe("a 205 → 195 Å"); // cubic → a
  });

  it("collapses a constant lattice to a single value", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Lamellar", lattice_d: 60 }]),
      member([{ phase: "Lamellar", lattice_d: 60 }]),
    ];
    const vars = ["a", "b"];
    const reading = seriesReading(members, (m) => vars[m.display_order]!);
    expect(reading.phases[0]!.latticeTrend).toBe("d = 60 Å"); // lamellar → d
  });

  it("emits coexistence + form-factor lines and omits no-phase members from indexed phases", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Pn3m", lattice_d: 195 }, { phase: "Lamellar", lattice_d: 60 }]),
      member([], "form_factor"),
    ];
    const vars = ["1:0.5", "1:1.5"];
    const reading = seriesReading(members, (m) => vars[m.display_order]!);
    expect(reading.coexistenceAt).toEqual(["1:0.5"]);
    expect(reading.formFactorOnlyAt).toEqual(["1:1.5"]);
    expect(new Set(reading.phases.map((p) => p.phase))).toEqual(new Set(["Pn3m", "Lamellar"]));
  });

  it("distinguishes a null member from a form-factor member (neither is an indexed phase)", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Pn3m", lattice_d: 200 }]),
      member([], "null"),
    ];
    const vars = ["x", "y"];
    const reading = seriesReading(members, (m) => vars[m.display_order]!);
    // A null member contributes no indexed phase and no form-factor line.
    expect(reading.formFactorOnlyAt).toEqual([]);
    expect(reading.phases.map((p) => p.phase)).toEqual(["Pn3m"]);
  });
});
