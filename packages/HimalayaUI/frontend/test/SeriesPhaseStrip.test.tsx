/**
 * SeriesPhaseStrip (Plan E, Task E-5) — derives PhaseSegment[] from member
 * snapshots (in variable / display order) and composes ui/PhaseStrip: one cell
 * per sample, coexistence gradient, hollow-dashed form-factor cell, distinct
 * null cell.
 */
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import type { SeriesMember, MemberSnapshotPhase, AssignmentState } from "../src/api";
import { SeriesPhaseStrip, segmentsFromMembers } from "../src/components/SeriesPhaseStrip";

let nextId = 1;
function member(phases: MemberSnapshotPhase[], state?: AssignmentState): SeriesMember {
  const id = nextId++;
  const st = state ?? (phases.length > 0 ? "indexed" : "form_factor");
  const primary = phases[0] ?? null;
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: {
      effective_peaks: [],
      confirmed_index: primary
        ? { id: id * 1000, phase: primary.phase, lattice_d: primary.lattice_d ?? 0,
            r_squared: 0.99, ngc: -1.5, peak_ids: [] }
        : null,
      confirmed_phases: phases,
      assignment_state: st,
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: null, created_at: null,
  };
}

describe("segmentsFromMembers", () => {
  it("maps N members to N segments in member order", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Pn3m", lattice_d: 205 }]),
      member([{ phase: "Pn3m", lattice_d: 195 }, { phase: "Lamellar", lattice_d: 60 }]),
      member([], "form_factor"),
      member([], "null"),
    ];
    const segs = segmentsFromMembers(members);
    expect(segs).toHaveLength(4);
    expect(segs[0]).toEqual({ phase: "Pn3m", coexistWith: null });
    expect(segs[1]!.coexistWith).toBe("Lamellar"); // coexistence gradient
    expect(segs[2]!.state).toBe("form_factor");
    expect(segs[3]!.state).toBe("null");
  });
});

describe("SeriesPhaseStrip", () => {
  it("renders one cell per member with the special-state cells", () => {
    nextId = 1;
    const members = [
      member([{ phase: "Pn3m", lattice_d: 205 }]),
      member([{ phase: "Pn3m", lattice_d: 195 }, { phase: "Lamellar", lattice_d: 60 }]),
      member([], "form_factor"),
      member([], "null"),
    ];
    render(<SeriesPhaseStrip members={members} />);
    const cells = screen.getAllByTestId("ps-seg");
    expect(cells).toHaveLength(4);
    expect(cells[2]).toHaveAttribute("data-state", "form_factor");
    expect(cells[3]).toHaveAttribute("data-state", "null");
  });
});
