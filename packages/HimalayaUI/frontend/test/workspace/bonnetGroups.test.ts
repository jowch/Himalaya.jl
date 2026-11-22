import { describe, expect, it } from "vitest";
import type { IndexEntry } from "../../src/api";
import { bonnetGroups, callRows } from "../../src/workspace/bonnetGroups";

let nextId = 1;
function ix(phase: string, a: number, peaks: readonly number[], score = 0.5): IndexEntry {
  return {
    id: nextId++, exposure_id: 1, phase, basis: 0.1, score, r_squared: null, lattice_d: a,
    ngc: null, status: "candidate", kind: "auto", inputs_hash: null, predicted_q: [],
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

  it("keeps lattices more than 1% off the Bonnet ratio apart", () => {
    expect(bonnetGroups([ix("Pn3m", 100, range(1, 4)), ix("Im3m", 127.9 * 1.011, range(5, 8))])).toEqual([]);
    expect(bonnetGroups([ix("Pn3m", 100, range(1, 4)), ix("Im3m", 127.9 * 1.009, range(5, 8))])).toHaveLength(1);
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
