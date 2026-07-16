import { describe, it, expect } from "vitest";
import type { SeriesMember, SeriesSample, CorpusSample } from "../../src/api";
import {
  membersToMemberData,
  recipeRowView,
  addableSamples,
  legendPhasesOf,
} from "../../src/print/pages/builderAdapters";

// ── fixtures ─────────────────────────────────────────────────────────────────

function member(
  id: number,
  snap: SeriesMember["snapshot"],
  overrides: Partial<SeriesMember> = {},
): SeriesMember {
  return {
    id,
    series_id: 10,
    exposure_id: id,
    display_order: id,
    band_height: 1,
    y_offset: 0,
    normalization: "none",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: snap,
    is_stale: false,
    created_by: null,
    created_at: null,
    ...overrides,
  };
}

function snap(
  phase: string | null,
  lattice_d: number | null = null,
  q1: number | null = null,
): SeriesMember["snapshot"] {
  if (phase === null) return null;
  const confirmedIndex =
    lattice_d !== null && q1 !== null
      ? {
          id: 1,
          phase,
          lattice_d,
          r_squared: 0.99,
          ngc: 3,
          peak_ids: [100],
        }
      : null;
  const effectivePeaks =
    q1 !== null
      ? [{ id: 100, q: q1, intensity: 1000, sharpness: 0.5, source: "auto" as const }]
      : [];
  return {
    effective_peaks: effectivePeaks,
    confirmed_index: confirmedIndex,
    confirmed_phases: confirmedIndex ? [{ phase, lattice_d }] : [],
    analysis_inputs_hash: "h",
  };
}

function seriesSample(id: number, sampleId: number, position: number): SeriesSample {
  return {
    id,
    series_id: 10,
    sample_id: sampleId,
    position,
    pinned: false,
    excluded: false,
  };
}

function corpusSample(id: number, name: string): CorpusSample {
  return {
    id,
    experiment_id: 1,
    name,
    notes: null,
    tags: [],
    q_units: "Å⁻¹",
  };
}

// ── membersToMemberData ───────────────────────────────────────────────────────

describe("membersToMemberData", () => {
  it("key is the member id as string", () => {
    const m = member(42, snap("Ia3d", 195, 0.057));
    const [row] = membersToMemberData([m]);
    expect(row!.key).toBe("42");
  });

  it("phases: single dominant phase when only confirmed_index", () => {
    const m = member(1, snap("Im3m", 150, 0.06));
    const [row] = membersToMemberData([m]);
    expect(row!.phases).toEqual(["Im3m"]);
  });

  it("phases: coexistence — confirmed_phases drives the list (dominant + coexist)", () => {
    const m = member(1, {
      effective_peaks: [{ id: 100, q: 0.057, intensity: 1000, sharpness: 0.5, source: "auto" }],
      confirmed_index: { id: 1, phase: "Ia3d", lattice_d: 195, r_squared: 0.99, ngc: 3, peak_ids: [100] },
      confirmed_phases: [
        { phase: "Ia3d", lattice_d: 195 },
        { phase: "Im3m", lattice_d: 150 },
      ],
      analysis_inputs_hash: "h",
    });
    const [row] = membersToMemberData([m]);
    expect(row!.phases).toEqual(["Ia3d", "Im3m"]);
  });

  it("phases: form_factor state → empty array", () => {
    const m = member(1, {
      effective_peaks: [],
      confirmed_index: null,
      assignment_state: "form_factor",
      analysis_inputs_hash: "h",
    });
    const [row] = membersToMemberData([m]);
    expect(row!.phases).toEqual([]);
  });

  it("phases: null snapshot → empty array", () => {
    const m = member(1, null);
    const [row] = membersToMemberData([m]);
    expect(row!.phases).toEqual([]);
  });

  it("variableValue: uses label_override when present", () => {
    const m = member(1, snap("Ia3d", 195, 0.057), { label_override: "1:0.5" });
    const [row] = membersToMemberData([m]);
    expect(row!.variableValue).toBe("1:0.5");
  });

  it("variableValue: falls back to exposure id label when no override", () => {
    // No label_override; exposure_id = 7
    const m = member(7, snap("Im3m", 150, 0.06));
    const [row] = membersToMemberData([m]);
    // Should produce the fallback label (matches toWaterfallRows label derivation)
    expect(typeof row!.variableValue).toBe("string");
    expect(String(row!.variableValue)).toContain("7");
  });

  it("dataLine: contains formatted lattice and q₁ values when available", () => {
    const m = member(1, snap("Ia3d", 195, 0.057));
    const [row] = membersToMemberData([m]);
    const line = String(row!.dataLine);
    expect(line).toContain("195");
    expect(line).toContain("0.057");
  });

  it("dataLine: empty string when no confirmed_index", () => {
    const m = member(1, null);
    const [row] = membersToMemberData([m]);
    expect(row!.dataLine).toBe("");
  });

  it("maps multiple members preserving order", () => {
    const members = [
      member(1, snap("Ia3d", 195, 0.057)),
      member(2, snap("Im3m", 150, 0.06)),
      member(3, null),
    ];
    const rows = membersToMemberData(members);
    expect(rows).toHaveLength(3);
    expect(rows[0]!.key).toBe("1");
    expect(rows[1]!.key).toBe("2");
    expect(rows[2]!.key).toBe("3");
  });
});

// ── recipeRowView ─────────────────────────────────────────────────────────────

describe("recipeRowView", () => {
  it("derives name from sampleNameById lookup", () => {
    const row = seriesSample(1, 5, 1);
    const nameById: Record<number, string> = { 5: "JC001-LipidA" };
    const view = recipeRowView(row, nameById);
    expect(view.name).toBe("JC001-LipidA");
  });

  it("falls back to 'Sample <id>' when the sample is not in the map", () => {
    const row = seriesSample(1, 99, 2);
    const view = recipeRowView(row, {});
    expect(view.name).toBe("Sample 99");
  });

  it("echoes the position field", () => {
    const row = seriesSample(1, 5, 3);
    const view = recipeRowView(row, { 5: "A" });
    expect(view.position).toBe(3);
  });
});

// ── addableSamples ────────────────────────────────────────────────────────────

describe("addableSamples", () => {
  it("excludes samples whose id is already in the draft recipe", () => {
    const corpus = [corpusSample(1, "A"), corpusSample(2, "B"), corpusSample(3, "C")];
    const draft: SeriesSample[] = [seriesSample(10, 1, 1), seriesSample(11, 3, 2)];
    const result = addableSamples(corpus, draft);
    expect(result.map((s) => s.id)).toEqual([2]);
  });

  it("returns all corpus samples when the draft is empty", () => {
    const corpus = [corpusSample(1, "A"), corpusSample(2, "B")];
    expect(addableSamples(corpus, [])).toHaveLength(2);
  });

  it("returns empty when all corpus samples are in the draft", () => {
    const corpus = [corpusSample(1, "A")];
    const draft: SeriesSample[] = [seriesSample(10, 1, 1)];
    expect(addableSamples(corpus, draft)).toHaveLength(0);
  });
});

// ── legendPhasesOf ────────────────────────────────────────────────────────────

describe("legendPhasesOf", () => {
  it("returns distinct non-null phases from the dominant phase of each row", () => {
    const members = [
      member(1, snap("Ia3d", 195, 0.057)),
      member(2, snap("Im3m", 150, 0.06)),
      member(3, snap("Ia3d", 195, 0.057)),  // duplicate
      member(4, null),                        // null phase
    ];
    const phases = legendPhasesOf(members);
    expect(phases).toHaveLength(2);
    expect(phases).toContain("Ia3d");
    expect(phases).toContain("Im3m");
  });

  it("includes coexistence phases from confirmed_phases", () => {
    const m = member(1, {
      effective_peaks: [{ id: 100, q: 0.057, intensity: 1000, sharpness: 0.5, source: "auto" }],
      confirmed_index: { id: 1, phase: "Ia3d", lattice_d: 195, r_squared: 0.99, ngc: 3, peak_ids: [100] },
      confirmed_phases: [
        { phase: "Ia3d", lattice_d: 195 },
        { phase: "Im3m", lattice_d: 150 },
      ],
      analysis_inputs_hash: "h",
    });
    const phases = legendPhasesOf([m]);
    expect(phases).toContain("Ia3d");
    expect(phases).toContain("Im3m");
  });

  it("returns empty array for all-null members", () => {
    expect(legendPhasesOf([member(1, null), member(2, null)])).toHaveLength(0);
  });

  it("returns empty array for empty input", () => {
    expect(legendPhasesOf([])).toHaveLength(0);
  });
});
