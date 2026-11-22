import { describe, expect, it } from "vitest";
import type { IndexEntry } from "../../src/api";
import { bonnetGroups, callRows, keyedRows } from "../../src/workspace/bonnetGroups";

let nextId = 1;
function ix(phase: string, a: number, peaks: readonly number[], score = 0.5): IndexEntry {
  return {
    id: nextId++, exposure_id: 1, phase, basis: 0.1, score, r_squared: null, lattice_d: a,
    status: "candidate", kind: "auto", inputs_hash: null, predicted_q: [],
    peaks: peaks.map((p, i) => ({ peak_id: p, ratio_position: i + 1, residual: 0, q_observed: p / 100 })),
  };
}
const range = (a: number, b: number) => Array.from({ length: b - a + 1 }, (_, i) => a + i);
const ids = (groups: readonly (readonly IndexEntry[])[]) => groups.map((g) => g.map((m) => m.id));

describe("bonnetGroups", () => {
  it("groups a Bonnet pair that shares peaks", () => {
    const c = [ix("Im3m", 127.9, range(1, 6), 0.9), ix("Pn3m", 100, [5, 6, 7, 8], 0.6)];
    expect(ids(bonnetGroups(c))).toEqual([[c[0]!.id, c[1]!.id]]);
  });

  it("does not group on a curvature match alone", () => {
    const c = [ix("Pn3m", 100, range(1, 6)), ix("Im3m", 127.9, [1, 2, 3, 7])];
    expect(bonnetGroups(c)).toEqual([]);
  });


  it("leaves out a third cubic that adds no peaks of its own", () => {
    const c = [ix("Pn3m", 100, range(1, 9), 0.8), ix("Im3m", 127.9, range(7, 16), 1), ix("Ia3d", 157.6, [2, 4, 9, 12], 0.9)];
    expect(ids(bonnetGroups(c))).toEqual([[c[0]!.id, c[1]!.id]]);
  });

  it("forms a triple when each cubic adds peaks", () => {
    const c = [ix("Pn3m", 100, range(1, 4)), ix("Im3m", 127.9, range(5, 8)), ix("Ia3d", 157.6, range(9, 12))];
    expect(ids(bonnetGroups(c))).toEqual([c.map((m) => m.id)]);
  });

  it("drops alias pairs that explain only peaks a better group already claims", () => {
    const c = [
      ix("Pn3m", 108, range(1, 9), 0.84), ix("Im3m", 138.1, range(7, 16), 1),
      ix("Pn3m", 76.4, [1, 3, 10, 12], 0.6), ix("Im3m", 97.7, [3, 5, 8, 14, 15], 0.78),
    ];
    expect(ids(bonnetGroups(c))).toEqual([[c[0]!.id, c[1]!.id]]);
  });

  it("drops a later group whose member explains only peaks an earlier group claimed", () => {
    const c = [
      ix("Pn3m", 100, range(1, 6), 0.9), ix("Im3m", 127.9, range(5, 10), 0.9),
      ix("Pn3m", 141.4, range(1, 4), 0.5), ix("Im3m", 180.9, [11, 12, 1, 2], 0.5),
    ];
    expect(ids(bonnetGroups(c))).toEqual([[c[0]!.id, c[1]!.id]]);
  });

  it("puts a candidate in at most one group", () => {
    const c = [ix("Pn3m", 100, range(1, 4)), ix("Im3m", 127.9, range(5, 10)), ix("Im3m", 128.5, range(11, 14))];
    expect(ids(bonnetGroups(c))).toEqual([[c[0]!.id, c[1]!.id]]);
  });

  it("never groups the same phase or a non-cubic phase", () => {
    expect(bonnetGroups([ix("Pn3m", 100, range(1, 4)), ix("Pn3m", 100, range(5, 8)), ix("Lamellar", 127.9, range(9, 12))])).toEqual([]);
  });

  it("ignores candidates without a lattice", () => {
    const c = [ix("Pn3m", 100, range(1, 4)), { ...ix("Im3m", 127.9, range(5, 8)), lattice_d: null }];
    expect(bonnetGroups(c)).toEqual([]);
  });

  it("counts a peak once whether it came from detection or a manual add", () => {
    const a = ix("Pn3m", 100, range(1, 4));
    const b = ix("Im3m", 127.9, range(3, 5));
    b.peaks = b.peaks.map((p) => ({ ...p, peak_id: p.peak_id + 1000 }));
    expect(bonnetGroups([a, b])).toEqual([]);
  });

  it("does not depend on the order of the candidates", () => {
    const c = [ix("Pn3m", 100, range(1, 4), 0.2), ix("Im3m", 127.9, range(5, 8), 0.9)];
    expect(ids(bonnetGroups([...c].reverse()))).toEqual(ids(bonnetGroups(c)));
  });
});

function fitted(phase: string, a: number, base: number, scatter: number, n = 4): IndexEntry {
  const predicted_q = Array.from({ length: n }, (_, i) => base * (i + 1));
  return {
    ...ix(phase, a, []),
    predicted_q,
    peaks: predicted_q.map((q, i) => ({
      peak_id: Math.round(base * 1e4) + i, ratio_position: i + 1, residual: 0,
      q_observed: q * (1 + (i % 2 === 0 ? scatter : -scatter)),
    })),
  };
}
const pair = (scatter: number, off: number) =>
  bonnetGroups([fitted("Pn3m", 100, 0.01, scatter), fitted("Im3m", 127.93 * (1 + off), 0.013, scatter)]);

describe("Bonnet tolerance from the lattice fits", () => {
  it("holds precise fits to the instrument floor", () => {
    expect(pair(0, 0.007)).toHaveLength(1);
    expect(pair(0, 0.011)).toEqual([]);
  });

  it("widens for scattered fits", () => {
    expect(pair(0.01, 0.015)).toHaveLength(1);
    expect(pair(0.01, 0.017)).toHaveLength(1);
    expect(pair(0.01, 0.02)).toEqual([]);
  });

  it("never exceeds the cap however scattered the fits", () => {
    expect(pair(0.05, 0.019)).toHaveLength(1);
    expect(pair(0.05, 0.021)).toEqual([]);
  });

  it("keeps out an alias that only the cap rejects", () => {
    const pn = fitted("Pn3m", 100, 0.01, 0.03);
    const im = fitted("Im3m", 100 * Math.sqrt(1.5), 0.013, 0.03);
    expect(bonnetGroups([pn, im])).toEqual([]);
    expect(bonnetGroups([pn, fitted("Im3m", 127.93, 0.013, 0.03)])).toHaveLength(1);
  });

  it("keeps out an alias inside any tolerance when it explains too few peaks of its own", () => {
    const c = [ix("Pn3m", 100, range(1, 6)), ix("Ia3d", 100 * Math.sqrt(2.5), [2, 4, 6, 7])];
    expect(bonnetGroups(c)).toEqual([]);
  });

  it("uses the floor alone for fits with fewer than three peaks", () => {
    const two = (phase: string, a: number, base: number) => fitted(phase, a, base, 0.05, 2);
    expect(bonnetGroups([two("Pn3m", 100, 0.01), two("Im3m", 127.93 * 1.007, 0.013)])).toHaveLength(1);
    expect(bonnetGroups([two("Pn3m", 100, 0.01), two("Im3m", 127.93 * 1.011, 0.013)])).toEqual([]);
  });
});

describe("callRows", () => {
  it("places a group where its best member ranks and lists members under it", () => {
    const lam = ix("Lamellar", 60, [20, 21]);
    const im = ix("Im3m", 127.9, range(5, 8));
    const hex = ix("Hexagonal", 70, [30, 31, 32]);
    const pn = ix("Pn3m", 100, range(1, 4));
    const rows = callRows([lam, im, hex, pn]);
    expect(rows.map((r) => (r.kind === "group" ? r.members.map((m) => m.id) : r.candidate.id))).toEqual([
      lam.id, [im.id, pn.id], hex.id,
    ]);
  });

  it("gives a group the number of peaks its members explain together", () => {
    const rows = callRows([ix("Im3m", 127.9, range(1, 6)), ix("Pn3m", 100, [5, 6, 7, 8])]);
    expect(rows[0]).toMatchObject({ kind: "group", peaks: 8 });
  });
});

describe("keyedRows", () => {
  it("numbers a group header before its members, in list order", () => {
    const lam = ix("Lamellar", 60, [20, 21]);
    const im = ix("Im3m", 127.9, range(5, 8));
    const pn = ix("Pn3m", 100, range(1, 4));
    const keyed = keyedRows(callRows([lam, im, pn]));
    expect(keyed.map((r) => (r.kind === "group" ? "group" : r.candidate.id))).toEqual([lam.id, "group", im.id, pn.id]);
  });
});
