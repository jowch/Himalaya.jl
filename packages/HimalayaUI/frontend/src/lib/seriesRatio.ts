/**
 * seriesRatio — format a phase's defining reflection ratio series as
 * "√2 : √3 : 2" (cubic) or "1 : 2 : 3" (lamellar), restricted to the
 * positions an index actually claims. Mirrors `src/phase.jl phaseratios`:
 * the per-phase array of radicands (h²+k²+l², the value UNDER the radical),
 * position-indexed by the 1-based `ratio_position` on each claimed peak.
 *
 * R4 finding L-9/L-10: the focus rail surfaces this ratio in the Phase-call
 * block (per active phase) and on each candidate row.
 */

// Radicands (values under √) per phase, 1-based by `ratio_position`.
// Source of truth: src/phase.jl `phaseratios` — keep these in the same order
// AND the same length, since `ratio_position` indexes straight into them.
// Hexagonal carries no 11: N = h²+hk+k² has no integer solution, so √11 is not
// a permitted 2D hexagonal reflection (#304).
const RADICANDS: Record<string, number[]> = {
  Lamellar:  [1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121],
  Hexagonal: [1, 3, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28],
  Square:    [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20],
  Pn3m:      [2, 3, 4, 6, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20, 21],
  Im3m:      [2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
  Ia3d:      [6, 8, 14, 16, 20, 22, 24, 26],
  Fm3m:      [3, 4, 8, 11, 12],
  Fd3m:      [3, 8, 11, 12, 16, 19, 24, 27, 32, 35, 36],
};

const MAX_TERMS = 4;

function term(radicand: number): string {
  const root = Math.sqrt(radicand);
  // Perfect square → bare integer (√4 → 2); otherwise the √N form.
  return Number.isInteger(root) ? String(root) : `√${radicand}`;
}

/** The single ratio term a peak at `position` (1-based `ratio_position`)
 *  carries for `phase`, e.g. "√2", "√3", "2". Used to label a candidate's
 *  claimed peaks on the trace (candidate hover). Empty when the phase or
 *  position is unknown. */
export function ratioTerm(phase: string, position: number): string {
  const radicands = RADICANDS[phase];
  if (!radicands || position < 1 || position > radicands.length) return "";
  return term(radicands[position - 1]!);
}

export function seriesRatio(phase: string, positions: number[]): string {
  const radicands = RADICANDS[phase];
  if (!radicands) return "";
  const claimed = [...new Set(positions)]
    .filter((p) => p >= 1 && p <= radicands.length)
    .sort((a, b) => a - b);
  if (claimed.length === 0) return "";
  const trimmed = claimed.slice(0, MAX_TERMS);
  const out = trimmed.map((p) => term(radicands[p - 1]!)).join(" : ");
  return claimed.length > MAX_TERMS ? `${out} …` : out;
}
