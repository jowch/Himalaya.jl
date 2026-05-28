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
import type { SeriesMember } from "../src/api";
import {
  colorFor,
  COMPARE_PALETTE,
  ORPHAN_FALLBACK,
  type GroupingMode,
} from "../src/lib/comparison/coloring";
import { phaseColor, PHASE_PALETTE } from "../src/phases";
import { angularHueDistance } from "../src/lib/color-distance";

function makeMember(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1,
    series_id: 100,
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
function sampleIdResolver(map: Record<number, number>): (m: SeriesMember) => number | null {
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

  it("COMPARE_PALETTE does not collide with the phase palette (string-equality)", () => {
    // First-line defence: no exact OKLCH string match. The hue-offset floor
    // below is the load-bearing perceptual check.
    const phaseColors = new Set(Object.values(PHASE_PALETTE));
    for (const c of COMPARE_PALETTE) {
      expect(phaseColors.has(c)).toBe(false);
    }
  });

  it("every COMPARE_PALETTE hue sits ≥13° from every PHASE_PALETTE hue", () => {
    // Perceptual floor — string-equality alone is the wrong check (round-1
    // review of #208/#251 caught a 3° hue collision that string-equality
    // missed). At L≈0.55 C≈0.13, 13° hue shift is ΔE2000 ≈ 2.5, the smallest
    // gap that keeps the byPhase/bySample modes from visually conflating.
    //
    // 13° (not 15°) is the floor because the eight phase hues pack the warm
    // and purple sectors tightly enough that ≥15° everywhere AND twelve
    // mutually-distinct entries is infeasible. Nine of the twelve sit at
    // exactly 15°; the 13° outlier is entry hue 5 vs Fm3m 18 (squeezed
    // between Fm3m and Hex). If the palette is re-laid, keep this floor
    // numeric, and update both the floor and the docstring in lockstep.
    const phaseHues = Object.values(PHASE_PALETTE).map(s => parseOklch(s).h);
    for (const c of COMPARE_PALETTE) {
      const ch = parseOklch(c).h;
      for (const ph of phaseHues) {
        const dist = angularHueDistance(ch, ph);
        expect(dist, `COMPARE ${c} (hue ${ch}) vs PHASE hue ${ph}`).toBeGreaterThanOrEqual(13);
      }
    }
  });
});

// --- Paper-tune verification (R8-N2 / round-2 finding) ----------------------
// Mirrors the methodology used in `phases.test.ts` for the R0b PHASE_PALETTE
// retune (#222 / #236): self-contained OKLCH→linear-sRGB→WCAG helper, then
// assert every palette entry clears AA (≥ 4.5:1) against `--plate`.
//
// The pre-R8-N2 palette was dark-tuned (L 0.76–0.80) for the retired dark
// surface; on the warm paper plate the bySample/distinct traces washed out.
// This block is what flips when the palette is retuned.

const PLATE = "oklch(0.992 0.004 90)";

function parseOklch(s: string): { L: number; C: number; h: number } {
  const m = /oklch\(\s*([\d.]+)\s+([\d.]+)\s+([\d.-]+)\s*\)/i.exec(s);
  if (!m) throw new Error(`not an oklch() string: ${s}`);
  return { L: parseFloat(m[1]!), C: parseFloat(m[2]!), h: parseFloat(m[3]!) };
}

function oklchToLinearSrgb({ L, C, h }: { L: number; C: number; h: number }): [number, number, number] {
  const hr = (h * Math.PI) / 180;
  const a = C * Math.cos(hr);
  const b = C * Math.sin(hr);
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.291485548 * b;
  const l = l_ * l_ * l_;
  const m = m_ * m_ * m_;
  const s = s_ * s_ * s_;
  return [
    +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s,
  ];
}

function relativeLuminance(s: string): number {
  const [r, g, b] = oklchToLinearSrgb(parseOklch(s)).map((c) => Math.min(1, Math.max(0, c)));
  return 0.2126 * r! + 0.7152 * g! + 0.0722 * b!;
}

function contrastRatio(a: string, b: string): number {
  const la = relativeLuminance(a);
  const lb = relativeLuminance(b);
  const lighter = Math.max(la, lb);
  const darker = Math.min(la, lb);
  return (lighter + 0.05) / (darker + 0.05);
}

describe("COMPARE_PALETTE — paper-tuned (R8-N2, round-2 finding)", () => {
  it("every palette colour passes WCAG AA (>= 4.5:1) on --plate", () => {
    for (const c of COMPARE_PALETTE) {
      const ratio = contrastRatio(c, PLATE);
      expect(ratio, `${c} on --plate`).toBeGreaterThanOrEqual(4.5);
    }
  });

  it("every palette colour sits in the paper-tune luminance band (0.50 ≤ L ≤ 0.58)", () => {
    // Same band the R0b PHASE_PALETTE retune (#222) sits in. A dark-tuned
    // entry (L 0.76–0.80) drifting in here is the exact regression R8-N2
    // calls out — fail loudly rather than wash out on paper.
    for (const c of COMPARE_PALETTE) {
      const { L } = parseOklch(c);
      expect(L, `${c} luminance out of paper-tune band`).toBeGreaterThanOrEqual(0.5);
      expect(L, `${c} luminance out of paper-tune band`).toBeLessThanOrEqual(0.58);
    }
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
