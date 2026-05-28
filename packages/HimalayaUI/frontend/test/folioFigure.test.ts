import { describe, it, expect } from "vitest";
import { buildPhaseStrip, buildMiniWaterfall } from "../src/lib/series/folioFigure";
import type { SeriesMember } from "../src/api";

function member(phase: string | null, peaks: number[] = [], order = 0): SeriesMember {
  return {
    id: order + 1,
    series_id: 1,
    exposure_id: order + 100,
    display_order: order,
    band_height: 1,
    y_offset: 0,
    normalization: "max",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot:
      phase === null && peaks.length === 0
        ? null
        : {
            effective_peaks: peaks.map((q, i) => ({
              id: i + 1,
              q,
              intensity: 1,
              sharpness: 1,
              source: "auto" as const,
            })),
            confirmed_index:
              phase === null
                ? null
                : {
                    id: order + 1,
                    phase,
                    lattice_d: 100,
                    r_squared: 0.99,
                    ngc: 3,
                    peak_ids: peaks.map((_, i) => i + 1),
                  },
            analysis_inputs_hash: "h",
          },
    is_stale: false,
    created_by: 1,
    created_at: null,
  };
}

describe("buildPhaseStrip", () => {
  it("emits one segment per member in order", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(s.segments).toHaveLength(2);
    expect(s.segments[0]!.phase).toBe("Pn3m");
    expect(s.segments[1]!.phase).toBe("Lamellar");
  });

  it("classifies a transition (first ≠ last)", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(s.kind).toBe("transition");
    expect(s.first).toBe("Pn3m");
    expect(s.last).toBe("Lamellar");
  });

  it("classifies a single phase throughout", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04], 0), member("Pn3m", [0.041], 1)]);
    expect(s.kind).toBe("throughout");
  });

  it("classifies no clear phase when no member is indexed", () => {
    const s = buildPhaseStrip([member(null, [], 0), member(null, [], 1)]);
    expect(s.kind).toBe("none");
    expect(s.segments[0]!.phase).toBeNull();
  });

  it("uses first/last INDEXED member for the caption, skipping unindexed ends", () => {
    const s = buildPhaseStrip([
      member(null, [], 0),
      member("Pn3m", [0.04], 1),
      member("Lamellar", [0.1], 2),
      member(null, [], 3),
    ]);
    expect(s.first).toBe("Pn3m");
    expect(s.last).toBe("Lamellar");
    expect(s.kind).toBe("transition");
  });

  it("renders a coexistence gradient when a resolver supplies a second phase", () => {
    const s = buildPhaseStrip([member("Pn3m", [0.04, 0.1], 0)], () => "Lamellar");
    expect(s.segments[0]!.coexistWith).toBe("Lamellar");
  });

  it("handles an empty member list", () => {
    const s = buildPhaseStrip([]);
    expect(s.segments).toHaveLength(0);
    expect(s.kind).toBe("none");
  });
});

describe("buildMiniWaterfall", () => {
  it("emits one trace path per member with a baseline y stacked top→bottom", () => {
    const wf = buildMiniWaterfall([member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]);
    expect(wf.traces).toHaveLength(2);
    // baselines descend down the SVG: later members sit lower (larger y)
    expect(wf.traces[1]!.baselineY).toBeGreaterThan(wf.traces[0]!.baselineY);
    expect(wf.traces[0]!.linePath.startsWith("M")).toBe(true);
    expect(wf.traces[0]!.color).toMatch(/oklch|var/);
  });

  it("places peak ticks at log-mapped x within the viewBox width", () => {
    const wf = buildMiniWaterfall([member("Pn3m", [0.04, 0.057], 0)]);
    const t = wf.traces[0]!;
    expect(t.ticks.length).toBeGreaterThan(0);
    for (const tk of t.ticks) {
      expect(tk.x).toBeGreaterThanOrEqual(0);
      expect(tk.x).toBeLessThanOrEqual(wf.width);
    }
  });

  it("scales viewBox height with member count", () => {
    const a = buildMiniWaterfall([member("Pn3m", [0.04], 0)]);
    const b = buildMiniWaterfall([
      member("Pn3m", [0.04], 0),
      member("Pn3m", [0.04], 1),
      member("Pn3m", [0.04], 2),
    ]);
    expect(b.height).toBeGreaterThan(a.height);
  });

  it("colors an unindexed trace with the neutral ink-faint token", () => {
    const wf = buildMiniWaterfall([member(null, [], 0)]);
    expect(wf.traces[0]!.color).toContain("ink-faint");
  });

  it("handles an empty member list without throwing", () => {
    const wf = buildMiniWaterfall([]);
    expect(wf.traces).toHaveLength(0);
    expect(wf.height).toBeGreaterThan(0);
  });
});
