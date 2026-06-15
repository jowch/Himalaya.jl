import { describe, it, expect } from "vitest";
import { buildSeriesFigureKeys } from "../../src/print/export/seriesFigureKeys";
import { COMPARE_PALETTE_HEX } from "../../src/print/export/traceColors";
import type { SeriesMember, MemberSnapshot, MemberSnapshotPhase } from "../../src/api";

function member(over: Partial<SeriesMember>, snap: Partial<MemberSnapshot> | null): SeriesMember {
  return {
    id: 1, series_id: 1, exposure_id: 100, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "max", color_override: null,
    label_override: null, q_window_min: null, q_window_max: null, peak_display: null,
    is_stale: false, created_by: null, created_at: null,
    snapshot: snap === null ? null : ({
      effective_peaks: [], confirmed_index: null, assignment_state: "indexed",
      analysis_inputs_hash: "h", ...snap,
    } as MemberSnapshot),
    ...over,
  };
}

const ctx = {
  sampleNameForExposure: (id: number | null) => (id === 100 ? "HEPES 37 °C" : id === 200 ? "HEPES 50 °C" : null),
  fallbackLabel: (m: SeriesMember) => `exp ${m.exposure_id}`,
  qUnits: "A-1",
};

const phase = (p: string, l: number | null): MemberSnapshotPhase => ({ phase: p, lattice_d: l });

describe("buildSeriesFigureKeys", () => {
  it("labels by sample name and walks the categorical palette by trace order", () => {
    const keys = buildSeriesFigureKeys(
      [
        member({ exposure_id: 100 }, { confirmed_phases: [phase("Im3m", 252)] }),
        member({ id: 2, exposure_id: 200, display_order: 1 }, { confirmed_phases: [phase("Pn3m", 197)] }),
      ],
      ctx,
    );
    expect(keys[0]!.label).toBe("HEPES 37 °C");
    expect(keys[1]!.label).toBe("HEPES 50 °C");
    expect(keys[0]!.color).toBe(COMPARE_PALETTE_HEX[0]);
    expect(keys[1]!.color).toBe(COMPARE_PALETTE_HEX[1]);
  });

  it("formats a cubic segment with lattice a and κ; a non-cubic gets d and no κ", () => {
    const keys = buildSeriesFigureKeys(
      [
        member({ exposure_id: 100 }, { confirmed_phases: [phase("Im3m", 252)] }),
        member({ id: 2, exposure_id: 200 }, { confirmed_phases: [phase("Lamellar", 65)] }),
      ],
      ctx,
    );
    expect(keys[0]!.segments[0]!.phase).toBe("Im3m");
    expect(keys[0]!.segments[0]!.detail).toMatch(/^a = 252 Å · κ = .+Å⁻²$/);
    expect(keys[1]!.segments[0]!.detail).toBe("d = 65 Å"); // no κ for lamellar
  });

  it("coexistence yields one segment per phase, equal billing (no dominant)", () => {
    const keys = buildSeriesFigureKeys(
      [member({ exposure_id: 100 }, { confirmed_phases: [phase("Im3m", 255), phase("Pn3m", 198)] })],
      ctx,
    );
    expect(keys[0]!.segments.map((s) => s.phase)).toEqual(["Im3m", "Pn3m"]);
    expect(keys[0]!.segments.every((s) => /κ/.test(s.detail))).toBe(true);
  });

  it("a form-factor member carries no segments and a phaseless note", () => {
    const keys = buildSeriesFigureKeys(
      [member({ exposure_id: 100 }, { assignment_state: "form_factor", confirmed_phases: [] })],
      ctx,
    );
    expect(keys[0]!.segments).toEqual([]);
    expect(keys[0]!.note).toBe("form factor (no Bragg)");
  });

  it("falls back to the member row label when the sample is unresolved", () => {
    const keys = buildSeriesFigureKeys(
      [member({ exposure_id: 999 }, { confirmed_phases: [phase("Im3m", 252)] })],
      ctx,
    );
    expect(keys[0]!.label).toBe("exp 999");
  });
});
