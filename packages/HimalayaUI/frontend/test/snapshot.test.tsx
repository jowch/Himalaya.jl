/**
 * Client-side `computeMemberSnapshot` (Plan §Phase 3, Task 3.3).
 *
 * Pins parity with the server-side `compute_member_snapshot` in
 * `packages/HimalayaUI/src/comparisons.jl`:
 * - effective_peaks: read from `queryKeys.peaks(exposureId)` cache, which
 *   already mirrors the server's `auto − exclude ∪ add` union shape per
 *   `routes_peaks.jl::GET /api/exposures/:id/peaks`.
 *   - Auto peaks: source="auto", `intensity` non-null, `sharpness` non-null.
 *   - Excluded auto peaks (`excluded=true`) are filtered out client-side.
 *   - Manual peaks: source="manual", `intensity: null`, `sharpness` is the
 *     sampled value from the trace (server fills it on insert).
 * - confirmed_index: highest-scored R²-passing index from the active custom
 *   group. R²-gate is 0.98 (mirrors Phase 2 backend constant
 *   `CONFIRMED_INDEX_R2_GATE`). Below 0.98 → null. Tie-break: score DESC,
 *   id ASC.
 * - analysis_inputs_hash: from the exposure entity cache.
 *
 * The single source of truth for the test fixtures' shape is `MemberSnapshot`
 * in `api.ts` — both the HTTP parser and the SSE handler land cache rows in
 * this shape.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { computeMemberSnapshot } from "../src/lib/comparison/snapshot";
import { queryKeys } from "../src/queries";
import type { Peak, IndexEntry, GroupEntry, Exposure } from "../src/api";

const EXPOSURE: Exposure = {
  id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
  status: null, image_path: null, image_version: "",
  trace_hash: null, analysis_inputs_hash: "sha256:abc",
  tags: [], sources: [],
};

const PEAK_AUTO_OK: Peak = {
  id: 10, exposure_id: 5, q: 0.123, intensity: 100.5, prominence: 0.5,
  sharpness: 30.0, source: "auto", excluded: false,
};
const PEAK_AUTO_EXCLUDED: Peak = {
  id: 11, exposure_id: 5, q: 0.456, intensity: 80.0, prominence: 0.3,
  sharpness: 20.0, source: "auto", excluded: true,
};
const PEAK_MANUAL: Peak = {
  id: 12, exposure_id: 5, q: 0.789, intensity: null, prominence: null,
  sharpness: 12.5, source: "manual", excluded: false,
};

function indexEntry(id: number, r2: number, score: number, peakIds: number[]): IndexEntry {
  return {
    id, exposure_id: 5, phase: "Pn3m", basis: 0.123,
    score, r_squared: r2, lattice_d: 50.0, ngc: -2.0,
    status: "candidate", kind: "auto", inputs_hash: "h1",
    peaks: peakIds.map((peak_id, i) => ({
      peak_id, ratio_position: i + 1, residual: 0.0, q_observed: 0.0,
    })),
    predicted_q: [],
  };
}

describe("computeMemberSnapshot", () => {
  let qc: QueryClient;
  beforeEach(() => { qc = new QueryClient(); });

  describe("effective_peaks", () => {
    it("filters out excluded auto peaks", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5),
        [PEAK_AUTO_OK, PEAK_AUTO_EXCLUDED, PEAK_MANUAL]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.effective_peaks.map((p) => p.id)).toEqual([10, 12]);
      // Sorted by q (mirrors server's `ORDER BY q`)
      expect(snap.effective_peaks.map((p) => p.q)).toEqual([0.123, 0.789]);
    });

    it("manual peaks keep intensity:null", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [PEAK_MANUAL]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.effective_peaks).toEqual([{
        id: 12, q: 0.789, intensity: null, sharpness: 12.5, source: "manual",
      }]);
    });

    it("auto peaks carry intensity + sharpness", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [PEAK_AUTO_OK]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.effective_peaks[0]).toEqual({
        id: 10, q: 0.123, intensity: 100.5, sharpness: 30.0, source: "auto",
      });
    });

    it("sharpness defaults to 0 when null on the cache row", () => {
      // Defensive: the MemberSnapshotPeak.sharpness type is `number` (not
      // nullable), but cache rows can carry null sharpness (rare auto-peak
      // edge cases or legacy data). The snapshot must coerce to a number
      // so the type contract holds.
      const odd: Peak = { ...PEAK_AUTO_OK, id: 99, sharpness: null };
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [odd]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.effective_peaks[0]!.sharpness).toBe(0);
    });

    it("returns empty array when peaks cache is missing", () => {
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.effective_peaks).toEqual([]);
    });
  });

  describe("confirmed_index (R²-gate + active custom group selection)", () => {
    it("picks highest-scored R²-passing index from active custom group", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [PEAK_AUTO_OK, PEAK_MANUAL]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.99, 0.85, [10, 12]),
        indexEntry(21, 0.99, 0.95, [10]),  // higher score
        indexEntry(22, 0.97, 0.99, [10]),  // below R²-gate
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true,
          members: [20, 21, 22] },
      ]);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.confirmed_index).not.toBeNull();
      expect(snap.confirmed_index!.id).toBe(21);
      expect(snap.confirmed_index!.peak_ids).toEqual([10]);
    });

    it("returns null when no index passes R²-gate (0.98)", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [PEAK_AUTO_OK]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.95, 0.99, [10]),
        indexEntry(21, 0.97, 0.99, [10]),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true, members: [20, 21] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index).toBeNull();
    });

    it("R²-gate is inclusive at 0.98", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.98, 0.5, []),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true, members: [20] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index?.id).toBe(20);
    });

    it("returns null when there is no active custom group", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.99, 0.99, [10]),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "auto", active: true, members: [20] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index).toBeNull();
    });

    it("returns null when the only custom group is inactive", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.99, 0.99, [10]),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: false, members: [20] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index).toBeNull();
    });

    it("only considers indices that are members of the active custom group", () => {
      // Index 22 has the highest score AND passes R², but it's NOT a member
      // of the active custom group — so it's excluded.
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), [PEAK_AUTO_OK]);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(20, 0.99, 0.85, [10]),
        indexEntry(22, 0.99, 0.99, [10]),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true, members: [20] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index!.id).toBe(20);
    });

    it("tie-break: same score → lower id wins (mirrors server `id ASC`)", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        indexEntry(25, 0.99, 0.85, []),
        indexEntry(20, 0.99, 0.85, []),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true,
          members: [20, 25] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index!.id).toBe(20);
    });

    it("score=null sorts last (mirrors server `score DESC NULLS LAST`)", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const withNullScore: IndexEntry = { ...indexEntry(20, 0.99, 0, []), score: null };
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [
        withNullScore,
        indexEntry(21, 0.99, 0.5, []),
      ]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true,
          members: [20, 21] },
      ]);
      expect(computeMemberSnapshot(5, qc).confirmed_index!.id).toBe(21);
    });

    it("confirmed_index carries id, phase, lattice_d, r_squared, ngc, peak_ids", () => {
      qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
      qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
      const ix: IndexEntry = {
        id: 20, exposure_id: 5, phase: "Im3m", basis: 0.5,
        score: 0.95, r_squared: 0.99, lattice_d: 75.5, ngc: -1.2,
        status: "candidate", kind: "auto", inputs_hash: "h",
        peaks: [
          { peak_id: 10, ratio_position: 1, residual: 0, q_observed: 0 },
          { peak_id: 12, ratio_position: 2, residual: 0, q_observed: 0 },
        ],
        predicted_q: [],
      };
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(5), [ix]);
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
        { id: 1, exposure_id: 5, kind: "custom", active: true, members: [20] },
      ]);
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.confirmed_index).toEqual({
        id: 20, phase: "Im3m", lattice_d: 75.5, r_squared: 0.99,
        ngc: -1.2, peak_ids: [10, 12],
      });
    });
  });

  describe("analysis_inputs_hash", () => {
    it("reads from exposure cache", () => {
      qc.setQueryData<Exposure>(queryKeys.exposure(5),
        { ...EXPOSURE, analysis_inputs_hash: "sha256:current" });
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.analysis_inputs_hash).toBe("sha256:current");
    });

    it("falls back to empty string when exposure cache or hash is missing", () => {
      // No exposure cache row. The MemberSnapshot type requires a string;
      // empty-string is the canonical "unknown" for the spec's drift detection
      // (snapshot:"" vs current:"sha256:..." → flagged stale, same as
      // the server's `is_member_stale`).
      const snap = computeMemberSnapshot(5, qc);
      expect(snap.analysis_inputs_hash).toBe("");
    });
  });
});
