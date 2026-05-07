/**
 * Coloring library tests (Plan §Phase 9, Task 9.1).
 *
 * Resolution order:
 *   member.color_override → grouping-mode default → fallback gray for orphans.
 *
 * Modes:
 *   bySample  — same `sample_id` (resolved via `sampleIdFor(member)`) → same color
 *   byPhase   — same confirmed_index.phase → same phase color from phases.ts;
 *               unindexed → fallback gray
 *   distinct  — palette walked by display_order; cycles when N > palette.length
 */
import { describe, it, expect } from "vitest";
import type { ComparisonMember } from "../src/api";
import {
  colorFor,
  COMPARE_PALETTE,
  ORPHAN_FALLBACK,
  type GroupingMode,
} from "../src/lib/comparison/coloring";
import { phaseColor } from "../src/phases";

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1,
    comparison_id: 100,
    exposure_id: 1,
    display_order: 0,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: {
      effective_peaks: [],
      confirmed_index: {
        id: 1, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.5,
        peak_ids: [],
      },
      analysis_inputs_hash: "abc",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

/** Trivial id-based sample resolver for tests: maps exposure_id to sample_id. */
function sampleIdResolver(map: Record<number, number>): (m: ComparisonMember) => number | null {
  return (m) => (m.exposure_id !== null ? map[m.exposure_id] ?? null : null);
}

describe("colorFor — palette + fallback", () => {
  it("exports a palette of 10–12 colors", () => {
    expect(COMPARE_PALETTE.length).toBeGreaterThanOrEqual(10);
    expect(COMPARE_PALETTE.length).toBeLessThanOrEqual(12);
  });

  it("ORPHAN_FALLBACK is a non-empty string", () => {
    expect(typeof ORPHAN_FALLBACK).toBe("string");
    expect(ORPHAN_FALLBACK.length).toBeGreaterThan(0);
  });
});

describe("colorFor — color_override always wins", () => {
  it("override beats bySample default", () => {
    const m = makeMember({ exposure_id: 1, color_override: "#abcdef" });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({ 1: 100 }),
    };
    expect(colorFor(m, "bySample", COMPARE_PALETTE, ctx)).toBe("#abcdef");
  });

  it("override beats byPhase default", () => {
    const m = makeMember({ color_override: "#deadbe" });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "byPhase", COMPARE_PALETTE, ctx)).toBe("#deadbe");
  });

  it("override beats distinct default", () => {
    const m = makeMember({ color_override: "#beefed" });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "distinct", COMPARE_PALETTE, ctx)).toBe("#beefed");
  });
});

describe("colorFor — bySample", () => {
  it("two members with the same sample_id get the same color", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1 });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({ 10: 7, 11: 7 }),
    };
    expect(colorFor(a, "bySample", COMPARE_PALETTE, ctx))
      .toBe(colorFor(b, "bySample", COMPARE_PALETTE, ctx));
  });

  it("two members with different sample_ids get different colors", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1 });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({ 10: 7, 11: 8 }),
    };
    expect(colorFor(a, "bySample", COMPARE_PALETTE, ctx))
      .not.toBe(colorFor(b, "bySample", COMPARE_PALETTE, ctx));
  });

  it("orphan exposure (NULL) → ORPHAN_FALLBACK", () => {
    const m = makeMember({ exposure_id: null });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "bySample", COMPARE_PALETTE, ctx)).toBe(ORPHAN_FALLBACK);
  });

  it("unresolved sample_id (cache miss) → ORPHAN_FALLBACK", () => {
    const m = makeMember({ exposure_id: 99 });
    const ctx = {
      allMembers: [m],
      // Resolver returns null for an unknown exposure id.
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "bySample", COMPARE_PALETTE, ctx)).toBe(ORPHAN_FALLBACK);
  });

  it("reordering doesn't shift colors (color is sample-property, not slot)", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1 });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({ 10: 7, 11: 8 }),
    };
    const beforeA = colorFor(a, "bySample", COMPARE_PALETTE, ctx);
    const reordered = {
      allMembers: [b, a],
      sampleIdFor: sampleIdResolver({ 10: 7, 11: 8 }),
    };
    const afterA = colorFor(a, "bySample", COMPARE_PALETTE, reordered);
    expect(afterA).toBe(beforeA);
  });
});

describe("colorFor — byPhase", () => {
  it("two members with the same phase get the same color (phaseColor)", () => {
    const a = makeMember({ id: 1, exposure_id: 1 });
    const b = makeMember({ id: 2, exposure_id: 2 });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(a, "byPhase", COMPARE_PALETTE, ctx)).toBe(phaseColor("Pn3m"));
    expect(colorFor(a, "byPhase", COMPARE_PALETTE, ctx))
      .toBe(colorFor(b, "byPhase", COMPARE_PALETTE, ctx));
  });

  it("members with different phases get different colors", () => {
    const a = makeMember({
      id: 1,
      snapshot: {
        effective_peaks: [],
        confirmed_index: { id: 1, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.5, peak_ids: [] },
        analysis_inputs_hash: "abc",
      },
    });
    const b = makeMember({
      id: 2,
      snapshot: {
        effective_peaks: [],
        confirmed_index: { id: 2, phase: "Hexagonal", lattice_d: 30, r_squared: 0.99, ngc: 0, peak_ids: [] },
        analysis_inputs_hash: "abc",
      },
    });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(a, "byPhase", COMPARE_PALETTE, ctx))
      .not.toBe(colorFor(b, "byPhase", COMPARE_PALETTE, ctx));
  });

  it("unindexed member (confirmed_index === null) → ORPHAN_FALLBACK", () => {
    const m = makeMember({
      snapshot: {
        effective_peaks: [],
        confirmed_index: null,
        analysis_inputs_hash: "abc",
      },
    });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "byPhase", COMPARE_PALETTE, ctx)).toBe(ORPHAN_FALLBACK);
  });

  it("snapshot null → ORPHAN_FALLBACK", () => {
    const m = makeMember({ snapshot: null });
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(m, "byPhase", COMPARE_PALETTE, ctx)).toBe(ORPHAN_FALLBACK);
  });
});

describe("colorFor — distinct", () => {
  it("each member gets a different palette color (in order)", () => {
    const members = Array.from({ length: 4 }, (_, i) =>
      makeMember({ id: i + 1, exposure_id: i + 1, display_order: i }),
    );
    const ctx = {
      allMembers: members,
      sampleIdFor: sampleIdResolver({}),
    };
    const colors = members.map((m) => colorFor(m, "distinct", COMPARE_PALETTE, ctx));
    expect(new Set(colors).size).toBe(4);
    expect(colors[0]).toBe(COMPARE_PALETTE[0]);
    expect(colors[1]).toBe(COMPARE_PALETTE[1]);
    expect(colors[2]).toBe(COMPARE_PALETTE[2]);
    expect(colors[3]).toBe(COMPARE_PALETTE[3]);
  });

  it("cycles palette when N > palette length", () => {
    const N = COMPARE_PALETTE.length + 2;
    const members = Array.from({ length: N }, (_, i) =>
      makeMember({ id: i + 1, exposure_id: i + 1, display_order: i }),
    );
    const ctx = {
      allMembers: members,
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(members[0]!, "distinct", COMPARE_PALETTE, ctx))
      .toBe(colorFor(members[COMPARE_PALETTE.length]!, "distinct", COMPARE_PALETTE, ctx));
    expect(colorFor(members[1]!, "distinct", COMPARE_PALETTE, ctx))
      .toBe(colorFor(members[COMPARE_PALETTE.length + 1]!, "distinct", COMPARE_PALETTE, ctx));
  });

  it("non-contiguous display_orders still cycle from index 0", () => {
    const a = makeMember({ id: 1, exposure_id: 1, display_order: 5 });
    const b = makeMember({ id: 2, exposure_id: 2, display_order: 17 });
    const ctx = {
      allMembers: [a, b],
      sampleIdFor: sampleIdResolver({}),
    };
    expect(colorFor(a, "distinct", COMPARE_PALETTE, ctx)).toBe(COMPARE_PALETTE[0]);
    expect(colorFor(b, "distinct", COMPARE_PALETTE, ctx)).toBe(COMPARE_PALETTE[1]);
  });
});

describe("colorFor — invalid mode falls through to ORPHAN_FALLBACK", () => {
  it("unknown mode → ORPHAN_FALLBACK", () => {
    const m = makeMember();
    const ctx = {
      allMembers: [m],
      sampleIdFor: sampleIdResolver({}),
    };
    // Force-cast an unknown mode to assert the fallback path.
    expect(colorFor(m, "garbage" as GroupingMode, COMPARE_PALETTE, ctx))
      .toBe(ORPHAN_FALLBACK);
  });
});
