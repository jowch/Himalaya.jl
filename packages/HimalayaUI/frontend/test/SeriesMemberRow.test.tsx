/**
 * SeriesMemberRow (Plan E, Task E-6) — per-member rail row: lattice (a/d, both
 * under coexistence), first-peak q₁ = min(effective_peaks.q), phase-coloured
 * names so coexistence rows self-decode; form-factor/null → "no Bragg peaks ·
 * q₁ —".
 */
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import type { SeriesMember, MemberSnapshotPhase, AssignmentState } from "../src/api";
import { SeriesMemberRow } from "../src/components/SeriesMemberRow";

let nextId = 1;
function member(
  phases: MemberSnapshotPhase[],
  observedQs: number[],
  state?: AssignmentState,
): SeriesMember {
  const id = nextId++;
  const st = state ?? (phases.length > 0 ? "indexed" : "form_factor");
  const primary = phases[0] ?? null;
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: {
      effective_peaks: observedQs.map((q, i) => ({
        id: id * 100 + i, q, intensity: 1, sharpness: 1, source: "auto" as const,
      })),
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

describe("SeriesMemberRow", () => {
  it("shows lattice a + first-peak q₁ for an indexed cubic member", () => {
    nextId = 1;
    const m = member([{ phase: "Pn3m", lattice_d: 205 }], [0.043, 0.061, 0.085]);
    render(<SeriesMemberRow member={m} variableLabel="1:0" />);
    const data = screen.getByTestId("member-row-data").textContent ?? "";
    expect(data).toContain("a = 205");
    expect(data).toContain("0.043"); // q₁ = min observed
    expect(screen.getByTestId("member-row-name")).toHaveTextContent("Pn3m");
  });

  it("shows both lattices for a coexistence member with both phase names", () => {
    nextId = 1;
    const m = member(
      [{ phase: "Pn3m", lattice_d: 195 }, { phase: "Lamellar", lattice_d: 60 }],
      [0.045, 0.105],
    );
    render(<SeriesMemberRow member={m} variableLabel="1:0.5" />);
    const data = screen.getByTestId("member-row-data").textContent ?? "";
    expect(data).toContain("a 195");
    expect(data).toContain("d 60");
    const names = screen.getByTestId("member-row-name").textContent ?? "";
    expect(names).toContain("Pn3m");
    expect(names).toContain("Lamellar");
  });

  it("renders the no-Bragg line for a form-factor member", () => {
    nextId = 1;
    const m = member([], [], "form_factor");
    render(<SeriesMemberRow member={m} variableLabel="1:1.5" />);
    const data = screen.getByTestId("member-row-data").textContent ?? "";
    expect(data).toMatch(/no bragg peaks/i);
    expect(data).toContain("q₁ —");
    expect(screen.getByTestId("member-row-name")).toHaveTextContent(/form factor/i);
  });

  it("phase names carry the phase colour via inline style", () => {
    nextId = 1;
    const m = member([{ phase: "Pn3m", lattice_d: 205 }], [0.043]);
    render(<SeriesMemberRow member={m} variableLabel="1:0" />);
    const name = screen.getByTestId("member-row-name").querySelector("[style]");
    expect(name?.getAttribute("style")).toMatch(/color/);
  });
});
