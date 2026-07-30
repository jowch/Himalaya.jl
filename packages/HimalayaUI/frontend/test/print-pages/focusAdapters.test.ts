import { describe, it, expect } from "vitest";
import type { Peak, IndexEntry, IndexPeakRef, Trace } from "../../src/api";
import {
  SYMS,
  customRefls,
  landsOn,
  latticeForFirstOrderOnPeak,
} from "../../src/lib/customIndex";
import { phaseColor } from "../../src/phases";
import {
  toTraceModel,
  peakClickAction,
  toDetectorRings,
  toCombSeries,
  customIndexPreview,
  CUSTOM_SYMS,
  buildDetectorCalibration,
} from "../../src/print/pages/focusAdapters";
import type { Experiment } from "../../src/api";

// ── fixtures ────────────────────────────────────────────────────────────────

function peak(over: Partial<Peak>): Peak {
  return {
    id: 1, exposure_id: 1, q: 0.1, intensity: 100, prominence: null,
    sharpness: null, source: "auto", excluded: false, ...over,
  };
}

describe("peakClickAction", () => {
  const auto = peak({ id: 7, source: "auto", excluded: false });
  const autoExcluded = peak({ id: 8, source: "auto", excluded: true });
  const manual = peak({ id: 9, source: "manual", excluded: false });

  it("a click on an auto peak toggles it excluded (disable a live one)", () => {
    expect(peakClickAction([auto, manual], 7)).toEqual({ kind: "toggle-exclude", excluded: true });
  });

  it("a click on an already-excluded auto peak restores it", () => {
    expect(peakClickAction([autoExcluded], 8)).toEqual({ kind: "toggle-exclude", excluded: false });
  });

  it("a click on a manual peak removes it", () => {
    expect(peakClickAction([auto, manual], 9)).toEqual({ kind: "remove" });
  });

  it("returns null for an unknown id (stale mark mid-reconcile)", () => {
    expect(peakClickAction([auto, manual], 999)).toBeNull();
  });
});

function ref(over: Partial<IndexPeakRef>): IndexPeakRef {
  return { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.1, ...over };
}

function ix(over: Partial<IndexEntry>): IndexEntry {
  return {
    id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.1, score: 0.9,
    r_squared: 0.99, lattice_d: 197, ngc: null, status: "candidate",
    kind: "auto", inputs_hash: null, peaks: [], predicted_q: [], ...over,
  };
}

// ── toTraceModel ──────────────────────────────────────────────────────────────

describe("toTraceModel", () => {
  const trace: Trace = { q: [0.1, 0.2], I: [1, 2], sigma: [0, 0] };

  it("maps each Peak to a PlotPeak (id/q/source/excluded) and passes phase through", () => {
    const peaks = [
      peak({ id: 5, q: 0.11, intensity: 50, source: "auto", excluded: false }),
      peak({ id: 6, q: 0.22, intensity: null, source: "manual", excluded: true }),
    ];
    const tm = toTraceModel(trace, peaks, "Im3m");
    expect(tm.trace).toBe(trace);
    expect(tm.phase).toBe("Im3m");
    expect(tm.peaks).toHaveLength(2);
    expect(tm.peaks[0]).toMatchObject({ id: 5, q: 0.11, intensity: 50, source: "auto", excluded: false });
    expect(tm.peaks[1]).toMatchObject({ id: 6, q: 0.22, source: "manual", excluded: true });
  });

  it("passes a null phase through unchanged", () => {
    expect(toTraceModel(trace, [], null).phase).toBeNull();
  });

  it("anchors a null-intensity (curation) peak to the trace curve at its q (nearest sample)", () => {
    // trace = q:[0.1,0.2] I:[1,2]; a manually-added peak at q=0.1 has no backend
    // intensity → nearest sample is q=0.1 → I=1. Without this it would drop to
    // baselineI in PlotPeaks and render flat on the axis.
    const tm = toTraceModel(trace, [peak({ id: 1, q: 0.1, intensity: null })], null);
    expect(tm.peaks[0]).toMatchObject({ intensity: 1 });
  });

  it("omits intensity (does not set it to undefined) when there is no curve to anchor to", () => {
    const empty: Trace = { q: [], I: [], sigma: [] };
    const tm = toTraceModel(empty, [peak({ id: 1, intensity: null })], null);
    expect("intensity" in tm.peaks[0]!).toBe(false);
  });
});

// ── toDetectorRings (FocusDetectorPanel 52-77) ────────────────────────────────

describe("toDetectorRings", () => {
  it("returns no rings and no phases when there are no active indices", () => {
    expect(toDetectorRings([], [peak({ id: 1, q: 0.1 })])).toEqual({
      rings: [],
      phases: [],
    });
  });

  it("claimed peaks -> coloured ring; predicted-absent -> ghost ring; leftover -> neutral", () => {
    const peaks = [
      peak({ id: 1, q: 0.10 }),
      peak({ id: 2, q: 0.20 }),
      peak({ id: 3, q: 0.50 }), // leftover (no index claims it)
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        // ratio_position joins each claim to its predicted_q slot (1-based);
        // position 3 (0.35) is unclaimed -> ghost.
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: 0.20 }),
        ],
        predicted_q: [0.10, 0.20, 0.35],
      }),
    ];
    const { rings, phases } = toDetectorRings(active, peaks);
    const c = phaseColor("Pn3m");
    // two claimed coloured rings
    expect(rings).toContainEqual({ q: 0.10, color: c });
    expect(rings).toContainEqual({ q: 0.20, color: c });
    // one ghost ring for the predicted-but-absent order
    expect(rings).toContainEqual({ q: 0.35, color: c, ghost: true });
    // leftover neutral ring is a bare q object (no color/ghost)
    expect(rings).toContainEqual({ q: 0.50 });
    expect(rings).toHaveLength(4);
    // the index emitted rings -> its phase is the caption identity
    expect(phases).toEqual(["Pn3m"]);
  });

  // FO-RESCORE2-P2 F2: "absent" means "predicted order whose ratio_position is
  // NOT claimed by this index" — the backend join — never a q-tolerance rescan
  // that can disagree with what the backend actually claimed.
  it("does NOT paint a ghost ring over a claimed order, even when its residual exceeds spanTol (both-rings repro)", () => {
    // peaks span [0.10, 0.30] → spanTol = 0.004. The index CLAIMS peak 2
    // (ratio_position 2, ideal 0.20) at q_observed 0.21 — 0.01 off, beyond tol.
    // The old tolerance scan saw no observed peak near 0.20 and painted BOTH an
    // observed ring at 0.21 AND a ghost at 0.20 for the same reflection.
    const peaks = [
      peak({ id: 1, q: 0.10 }),
      peak({ id: 2, q: 0.21 }),
      peak({ id: 3, q: 0.30 }), // leftover, widens the span
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: 0.21, residual: 0.01 }),
        ],
        predicted_q: [0.10, 0.20],
      }),
    ];
    const { rings } = toDetectorRings(active, peaks);
    const c = phaseColor("Pn3m");
    // the claimed reflection renders ONCE, at its observed q
    expect(rings).toContainEqual({ q: 0.21, color: c });
    // …and never also as a ghost at its predicted q
    expect(rings).not.toContainEqual({ q: 0.20, color: c, ghost: true });
    // no ghost at all: both predicted orders are claimed
    expect(rings.filter((r) => typeof r === "object" && r.ghost)).toHaveLength(0);
  });

  it("ghosts exactly the UNCLAIMED predicted orders of a claiming index", () => {
    const peaks = [peak({ id: 1, q: 0.10 }), peak({ id: 2, q: 0.30 })];
    const active = [
      ix({
        id: 1, phase: "Im3m",
        peaks: [ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 })],
        predicted_q: [0.10, 0.14, 0.17], // positions 2 and 3 unclaimed
      }),
    ];
    const { rings } = toDetectorRings(active, peaks);
    const c = phaseColor("Im3m");
    expect(rings).toContainEqual({ q: 0.14, color: c, ghost: true });
    expect(rings).toContainEqual({ q: 0.17, color: c, ghost: true });
    expect(rings).not.toContainEqual({ q: 0.10, color: c, ghost: true });
  });

  // FO-RING caption identity: `phases` names only the phases that actually put
  // a ring (coloured or ghost) on the frame — derived from the SAME walk that
  // emits the rings, never from bare ix.phase.
  describe("phases (ring-identity caption source)", () => {
    it("EXCLUDES a fully-landed custom index (peaks: [], every predicted_q matched -> zero rings)", () => {
      // insert_custom_index! writes no index_peaks rows, so a committed custom
      // index arrives with peaks: []. When its fit is perfect every predicted_q
      // sits within tol of an observed peak, no ghost survives, and it colours
      // ZERO rings (its peaks render as neutral leftovers). The caption must
      // not chip it.
      const peaks = [peak({ id: 1, q: 0.10 }), peak({ id: 2, q: 0.20 })];
      const active = [
        ix({
          id: 1, phase: "Im3m",
          peaks: [ref({ peak_id: 1, q_observed: 0.10 })],
          predicted_q: [0.10],
        }),
        ix({
          id: 2, phase: "Pn3m", kind: "speculative",
          peaks: [], // landed custom: pure lattice hypothesis, no claimed peaks
          predicted_q: [0.10, 0.20], // both within tol of observed -> no ghosts
        }),
      ];
      const { rings, phases } = toDetectorRings(active, peaks);
      // the sibling with claimed peaks IS the caption; the ring-less custom is not
      expect(phases).toEqual(["Im3m"]);
      // sanity: the custom emitted no Pn3m-coloured ring (solid or ghost)
      expect(
        rings.filter((r) => typeof r === "object" && r.color === phaseColor("Pn3m")),
      ).toHaveLength(0);
    });

    it("INCLUDES a peaks-empty index that emits at least one ghost (unmatched predicted_q)", () => {
      // A ghost ring still carries the phase hue, so it counts as identity on
      // the frame.
      const peaks = [peak({ id: 1, q: 0.10 }), peak({ id: 2, q: 0.20 })];
      const active = [
        ix({
          id: 1, phase: "Pn3m", kind: "speculative",
          peaks: [],
          predicted_q: [0.10, 0.35], // 0.35 unmatched -> one ghost ring
        }),
      ];
      const { rings, phases } = toDetectorRings(active, peaks);
      expect(phases).toEqual(["Pn3m"]);
      expect(rings).toContainEqual({ q: 0.35, color: phaseColor("Pn3m"), ghost: true });
    });

    it("dedupes a phase carried by two ring-emitting indices, keeping first position", () => {
      const peaks = [peak({ id: 1, q: 0.10 }), peak({ id: 2, q: 0.20 }), peak({ id: 3, q: 0.30 })];
      const active = [
        ix({ id: 1, phase: "Pn3m", peaks: [ref({ peak_id: 1, q_observed: 0.10 })], predicted_q: [0.10] }),
        ix({ id: 2, phase: "Im3m", peaks: [ref({ peak_id: 2, q_observed: 0.20 })], predicted_q: [0.20] }),
        ix({ id: 3, phase: "Pn3m", peaks: [ref({ peak_id: 3, q_observed: 0.30 })], predicted_q: [0.30] }),
      ];
      expect(toDetectorRings(active, peaks).phases).toEqual(["Pn3m", "Im3m"]);
    });
  });
});

// ── toCombSeries (CombPanel) ──────────────────────────────────────────────────

describe("toCombSeries", () => {
  it("builds a phase row with color, √N teeth labels, and a leftover set", () => {
    // Pn3m predicted q for a lattice that puts √2 on 0.10:
    const val = latticeForFirstOrderOnPeak("Pn3m", 0.10);
    const refls = customRefls("Pn3m", val); // N = [2,3,4,6,8,9]
    const q2 = refls[0]!.q; // ~0.10 (√2)
    const q3 = refls[1]!.q; // √3

    const peaks = [
      peak({ id: 1, q: q2 }),
      peak({ id: 2, q: q3 }),
      peak({ id: 9, q: 0.99 }), // leftover (no index claims it)
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m", lattice_d: 197, r_squared: 0.995,
        basis: q2,
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: q2, residual: 0.0 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: q3, residual: 0.0 }),
        ],
        predicted_q: [q2, q3],
      }),
    ];

    const { assigned, leftover } = toCombSeries(active, peaks);
    expect(assigned).toHaveLength(1);
    const s = assigned[0]!;
    expect(s.phase).toBe("Pn3m");
    expect(s.color).toBe(phaseColor("Pn3m"));
    expect(s.rSquared).toBe(0.995);
    expect(s.latticeLabel).toBe("a = 197 Å");
    // teeth labels are "√N" strings derived from the nearest reflection
    expect(s.teeth.map((t) => t.label)).toEqual(["√2", "√3"]);
    // both predicted teeth are claimed by an observed peak -> observed true
    expect(s.teeth.every((t) => t.observed)).toBe(true);
    // leftover = the unclaimed observed peak q
    expect(leftover).toEqual([0.99]);
  });

  it("marks a predicted-but-absent order observed:false (no claimed peak at that q)", () => {
    const val = latticeForFirstOrderOnPeak("Pn3m", 0.10);
    const refls = customRefls("Pn3m", val);
    const q2 = refls[0]!.q;
    const q3 = refls[1]!.q;
    const peaks = [peak({ id: 1, q: q2 })]; // only √2 observed
    const active = [
      ix({
        id: 1, phase: "Pn3m", lattice_d: null, r_squared: null, basis: q2,
        peaks: [ref({ peak_id: 1, ratio_position: 1, q_observed: q2, residual: 0 })],
        predicted_q: [q2, q3], // q3 predicted but absent
      }),
    ];
    const { assigned, leftover } = toCombSeries(active, peaks);
    const s = assigned[0]!;
    expect(s.teeth.find((t) => t.label === "√2")!.observed).toBe(true);
    expect(s.teeth.find((t) => t.label === "√3")!.observed).toBe(false);
    expect(s.latticeLabel).toBeUndefined();
    expect(s.rSquared).toBeUndefined();
    expect(leftover).toEqual([]);
  });

  // FO-RESCORE2-P2 F2: the cart counts ix.peaks.length "reflections"; the comb
  // must light the SAME reflections. A tooth is observed iff a claimed
  // IndexPeakRef joins to it by ratio_position (1-based index into
  // predicted_q) — never by re-matching q within a second tolerance regime.
  it("renders a claimed peak observed even when its q_observed is farther from the tooth than spanTol (F2 repro)", () => {
    // peaks span [0.10, 0.30] → spanTol = 0.004. The backend claimed peak 2
    // for ratio_position 2 (ideal 0.20) at q_observed 0.21 — 0.01 off, well
    // beyond tol. The old tolerance rescan rejected the match: the cart said
    // "2 reflections" while the comb showed 1 observed tooth.
    const peaks = [
      peak({ id: 1, q: 0.10 }),
      peak({ id: 2, q: 0.21 }),
      peak({ id: 3, q: 0.30 }), // leftover, widens the span
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: 0.21, residual: 0.01 }),
        ],
        predicted_q: [0.10, 0.20],
      }),
    ];
    const { assigned } = toCombSeries(active, peaks);
    const s = assigned[0]!;
    expect(s.teeth[1]!.observed).toBe(true);
    expect(s.teeth[1]!.residual).toBeCloseTo((0.21 - 0.20) / 0.20, 6); // 0.05
  });

  it("observed-teeth count equals ix.peaks.length (cart–comb consistency invariant)", () => {
    const peaks = [
      peak({ id: 1, q: 0.10 }),
      peak({ id: 2, q: 0.145 }),
      peak({ id: 3, q: 0.21 }),
      peak({ id: 4, q: 0.40 }), // unclaimed → leftover
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: 0.145, residual: 0.005 }),
          ref({ peak_id: 3, ratio_position: 4, q_observed: 0.21, residual: 0.01 }),
        ],
        predicted_q: [0.10, 0.14, 0.17, 0.20], // position 3 unclaimed
      }),
    ];
    const { assigned } = toCombSeries(active, peaks);
    const s = assigned[0]!;
    expect(s.teeth.filter((t) => t.observed)).toHaveLength(active[0]!.peaks.length);
    expect(s.teeth.map((t) => t.observed)).toEqual([true, true, false, true]);
  });

  it("derives the residual from the TRUE ideal (predicted_q[rpos-1]); a peak below prediction gets a negative fraction", () => {
    // Backend residual is abs(q_obs − ideal) (pipeline.jl), so reconstructing
    // the ideal as q_obs − residual is wrong whenever the peak sits BELOW the
    // prediction. q_observed = 0.198 vs ideal 0.20 → fraction −0.01.
    const peaks = [
      peak({ id: 1, q: 0.10 }),
      peak({ id: 2, q: 0.198 }),
      peak({ id: 3, q: 0.50 }), // widens span so old tol would have matched too
    ];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 }),
          ref({ peak_id: 2, ratio_position: 2, q_observed: 0.198, residual: 0.002 }),
        ],
        predicted_q: [0.10, 0.20],
      }),
    ];
    const { assigned } = toCombSeries(active, peaks);
    expect(assigned[0]!.teeth[1]!.observed).toBe(true);
    expect(assigned[0]!.teeth[1]!.residual).toBeCloseTo(-0.01, 6);
  });

  it("ignores a claimed ref whose ratio_position is out of bounds of predicted_q (defensive, no crash)", () => {
    const peaks = [peak({ id: 1, q: 0.10 })];
    const active = [
      ix({
        id: 1, phase: "Pn3m",
        peaks: [
          ref({ peak_id: 1, ratio_position: 1, q_observed: 0.10, residual: 0 }),
          ref({ peak_id: 99, ratio_position: 7, q_observed: 0.26, residual: 0 }), // out of range
        ],
        predicted_q: [0.10, 0.20],
      }),
    ];
    const { assigned } = toCombSeries(active, peaks);
    expect(assigned[0]!.teeth.map((t) => t.observed)).toEqual([true, false]);
  });

  it("does not crash on an unknown symmetry; labels fall back gracefully", () => {
    const active = [
      ix({ id: 1, phase: "Square", basis: 0.1, predicted_q: [0.1, 0.2], peaks: [] }),
    ];
    const { assigned } = toCombSeries(active, [peak({ id: 1, q: 0.1 })]);
    const s = assigned[0]!;
    expect(s.phase).toBe("Square");
    expect(s.teeth).toHaveLength(2);
    // fallback labels are non-empty strings, no throw
    expect(s.teeth.every((t) => typeof t.label === "string" && t.label.length > 0)).toBe(true);
  });

  it("labels a full Pn3m index BEYOND SYMS.Ms by the q-law (no √9 clamp on √10..√16)", () => {
    // q ∝ √N for cubic, anchored on √2 (basis = q0). Build predicted_q from a
    // single basis so q/q0 ratios are exact: q_N = q0·√(N/2).
    const q0 = 0.10; // √2 order
    const Ns = [2, 3, 4, 6, 8, 9, 10, 11, 12, 14, 16];
    const predicted_q = Ns.map((N) => q0 * Math.sqrt(N / 2));
    const active = [
      ix({ id: 1, phase: "Pn3m", basis: q0, lattice_d: null, r_squared: null,
           predicted_q, peaks: [] }),
    ];
    const { assigned } = toCombSeries(active, []);
    const labels = assigned[0]!.teeth.map((t) => t.label);
    expect(labels).toEqual(Ns.map((N) => `√${N}`));
    // the tail orders that USED to clamp to √9 are now correct:
    expect(labels.slice(-3)).toEqual(["√12", "√14", "√16"]);
  });

  it("labels Hexagonal orders beyond SYMS.Ms (no clamp to the last drawn order)", () => {
    // hex: q ∝ √M, anchored on M=1 (n1=1). q_M = q0·√M.
    const q0 = 0.05; // M=1
    // The backend series past SYMS.Hexagonal.Ms = [1,3,4,7,9,12]. Labels are
    // derived by inverting the q-law, not by indexing Ms, so the tail orders
    // must come out exact rather than clamping to √12.
    // (Pre-#304 this case used √11, which the core series wrongly carried at
    // position 6; it is not a permitted 2D hexagonal reflection.)
    const Ms = [1, 3, 4, 7, 9, 12, 13, 16];
    const predicted_q = Ms.map((M) => q0 * Math.sqrt(M));
    const active = [
      ix({ id: 1, phase: "Hexagonal", basis: q0, lattice_d: null, r_squared: null,
           predicted_q, peaks: [] }),
    ];
    const { assigned } = toCombSeries(active, []);
    const labels = assigned[0]!.teeth.map((t) => t.label);
    expect(labels).toEqual(["√1", "√3", "√4", "√7", "√9", "√12", "√13", "√16"]);
    expect(labels.slice(-2)).toEqual(["√13", "√16"]);
  });

  it("labels a Lamellar index by the linear q-law (q ∝ N) beyond SYMS.Ms", () => {
    // lamellar: q ∝ N, n1=1. orders 1..8 (SYMS.Ms stops at 5).
    const q0 = 0.04;
    const Ns = [1, 2, 3, 4, 5, 6, 7, 8];
    const predicted_q = Ns.map((N) => q0 * N);
    const active = [
      ix({ id: 1, phase: "Lamellar", basis: q0, lattice_d: null, r_squared: null,
           predicted_q, peaks: [] }),
    ];
    const { assigned } = toCombSeries(active, []);
    expect(assigned[0]!.teeth.map((t) => t.label)).toEqual(Ns.map((N) => `√${N}`));
  });
});

// ── CUSTOM_SYMS ───────────────────────────────────────────────────────────────

describe("CUSTOM_SYMS", () => {
  it("exposes one entry per modal symmetry with picker metadata derived from SYMS", () => {
    const names = CUSTOM_SYMS.map((s) => s.name);
    expect(names).toEqual([
      "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m", "Lamellar", "Hexagonal", "Square",
    ]);
    const pn3m = CUSTOM_SYMS.find((s) => s.name === "Pn3m")!;
    expect(pn3m.min).toBe(SYMS.Pn3m!.min);
    expect(pn3m.max).toBe(SYMS.Pn3m!.max);
    expect(pn3m.paramName).toBe("a");
    expect(pn3m.unit).toBe("Å");
    const lam = CUSTOM_SYMS.find((s) => s.name === "Lamellar")!;
    expect(lam.paramName).toBe("d");
  });
});

// ── customIndexPreview ────────────────────────────────────────────────────────

describe("customIndexPreview", () => {
  it("builds a preview series + a 'K of M reflections land' fit", () => {
    const sym = "Pn3m";
    const val = latticeForFirstOrderOnPeak(sym, 0.10); // first order on 0.10
    const refls = customRefls(sym, val); // N = [2,3,4,6,8,9]
    // observe the first two reflection q's exactly -> 2 land.
    const observed = [refls[0]!.q, refls[1]!.q];

    const { previewSeries, fit } = customIndexPreview(sym, val, observed);
    expect(previewSeries.phase).toBe(sym);
    expect(previewSeries.color).toBe(phaseColor(sym));
    // total = number of reflections
    expect(fit.total).toBe(refls.length);
    // landed matches the lib's landsOn count (= 2 here)
    expect(fit.landed).toBe(landsOn(sym, val, observed));
    expect(fit.landed).toBe(2);

    // teeth: labels √N, the two matched teeth observed:true with a residual,
    // the rest observed:false.
    expect(previewSeries.teeth).toHaveLength(refls.length);
    expect(previewSeries.teeth[0]!.label).toBe("√2");
    expect(previewSeries.teeth[0]!.observed).toBe(true);
    expect(previewSeries.teeth[0]!.residual).toBeCloseTo(0, 6);
    expect(previewSeries.teeth[1]!.observed).toBe(true);
    // an unobserved higher order
    expect(previewSeries.teeth[2]!.observed).toBe(false);
    expect(previewSeries.teeth[2]!.residual).toBeUndefined();
  });

  it("reports a non-zero fractional residual for a slightly-off observed peak", () => {
    const sym = "Pn3m";
    const val = latticeForFirstOrderOnPeak(sym, 0.10);
    const refls = customRefls(sym, val);
    const q1 = refls[0]!.q;
    const off = q1 * 1.01; // 1% high, inside the 0.022 tol
    const { previewSeries, fit } = customIndexPreview(sym, val, [off]);
    expect(fit.landed).toBe(1);
    const t0 = previewSeries.teeth[0]!;
    expect(t0.observed).toBe(true);
    expect(t0.residual).toBeCloseTo(0.01, 4);
  });

  it("returns empty teeth and a zero fit for an unknown symmetry", () => {
    const { previewSeries, fit } = customIndexPreview("Nonsense", 200, [0.1]);
    expect(previewSeries.teeth).toEqual([]);
    expect(fit).toEqual({ landed: 0, total: 0 });
  });
});

// ── buildDetectorCalibration ──────────────────────────────────────────────────

function makeExp(over: Partial<Experiment> = {}): Experiment {
  return {
    id: 1, name: "E", path: "/p", data_dir: "d", analysis_dir: "a",
    manifest_path: null, created_at: "", q_units: "A-1",
    beam_center_x: 420.791, beam_center_y: 838.83, pixel_size_um: 172,
    energy_kev: 12, flight_path_m: 2.5, ...over,
  };
}

describe("buildDetectorCalibration", () => {
  it("builds a full calibration with µm→mm and m→mm conversions", () => {
    const cal = buildDetectorCalibration(makeExp(), { w: 1475, h: 1679 });
    expect(cal).toEqual({
      beamCenterPx: { x: 420.791, y: 838.83 },
      imageSizePx: { w: 1475, h: 1679 },
      sampleDistanceMm: 2500,
      pixelSizeMm: 0.172,
      energyKeV: 12,
    });
  });

  it("returns null when rawSize is missing", () => {
    expect(buildDetectorCalibration(makeExp(), null)).toBeNull();
  });

  it("returns null when experiment is undefined", () => {
    expect(buildDetectorCalibration(undefined, { w: 100, h: 100 })).toBeNull();
  });

  it.each(["beam_center_x", "beam_center_y", "pixel_size_um", "energy_kev", "flight_path_m"] as const)(
    "returns null when %s is null", (field) => {
      expect(buildDetectorCalibration(makeExp({ [field]: null }), { w: 100, h: 100 })).toBeNull();
    });

  it("returns null for non-positive or non-finite dims", () => {
    expect(buildDetectorCalibration(makeExp(), { w: 0, h: 100 })).toBeNull();
    expect(buildDetectorCalibration(makeExp(), { w: NaN, h: 100 })).toBeNull();
  });
});
