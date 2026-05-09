/**
 * prefetchColdMembers — extracted helper that warms the four TanStack
 * cache keys (exposure / peaks / indices / groups) required by
 * computeMemberSnapshot for any member whose corresponding key is missing.
 * Per-key cold detection: only the missing keys are prefetched.
 *
 * Two call sites depend on this:
 *   - ComparePageEdit.handleSave (PR #49)
 *   - ConflictModal.buildOverwritePayload (PR #92, fix #74)
 * Both must move in lockstep when computeMemberSnapshot's required-key
 * set changes — that's why this helper exists. Issue #93.
 */
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { prefetchColdMembers } from "../src/lib/comparison/prefetchMembers";
import { queryKeys } from "../src/queries";
import * as api from "../src/api";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
    },
  });
}

let exposureSpy: ReturnType<typeof vi.spyOn>;
let peaksSpy: ReturnType<typeof vi.spyOn>;
let indicesSpy: ReturnType<typeof vi.spyOn>;
let groupsSpy: ReturnType<typeof vi.spyOn>;

beforeEach(() => {
  // Resolve to minimal valid shapes — the test asserts call counts, not
  // payload contents.
  exposureSpy = vi.spyOn(api, "getExposure").mockResolvedValue({
    id: 0, sample_id: 0, name: "x", file_path: "/x", selected: false,
    rejected: false, rejection_reason: null, exposure_type: "simple",
    notes: null, tags: [], analysis_inputs_hash: null, status: "ok",
  } as unknown as api.Exposure);
  peaksSpy   = vi.spyOn(api, "listPeaks").mockResolvedValue([]);
  indicesSpy = vi.spyOn(api, "listIndices").mockResolvedValue([]);
  groupsSpy  = vi.spyOn(api, "listGroups").mockResolvedValue([]);
});
afterEach(() => { vi.restoreAllMocks(); });

describe("prefetchColdMembers", () => {
  it("no-op when all four keys are warm for every member", async () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.exposure(100), { id: 100 });
    qc.setQueryData(queryKeys.peaks(100), []);
    qc.setQueryData(queryKeys.indices(100), []);
    qc.setQueryData(queryKeys.groups(100), []);

    await prefetchColdMembers([100], qc);

    expect(exposureSpy).not.toHaveBeenCalled();
    expect(peaksSpy).not.toHaveBeenCalled();
    expect(indicesSpy).not.toHaveBeenCalled();
    expect(groupsSpy).not.toHaveBeenCalled();
  });

  it("no-op when the input list is empty", async () => {
    const qc = makeQc();
    await prefetchColdMembers([], qc);
    expect(exposureSpy).not.toHaveBeenCalled();
    expect(peaksSpy).not.toHaveBeenCalled();
    expect(indicesSpy).not.toHaveBeenCalled();
    expect(groupsSpy).not.toHaveBeenCalled();
  });

  it("fires all four fetches for a fully-cold exposure", async () => {
    const qc = makeQc();

    await prefetchColdMembers([100], qc);

    expect(exposureSpy).toHaveBeenCalledWith(100);
    expect(peaksSpy).toHaveBeenCalledWith(100);
    expect(indicesSpy).toHaveBeenCalledWith(100);
    expect(groupsSpy).toHaveBeenCalledWith(100);
  });

  it("per-key cold detection — only the missing keys are fetched", async () => {
    // Issue #93 bonus: tighten cold detection from all-or-nothing per
    // exposure to per-key. If exposure / peaks / indices are warm but
    // groups is cold, only listGroups should fire.
    const qc = makeQc();
    qc.setQueryData(queryKeys.exposure(100), { id: 100 });
    qc.setQueryData(queryKeys.peaks(100), []);
    qc.setQueryData(queryKeys.indices(100), []);
    // groups intentionally not seeded.

    await prefetchColdMembers([100], qc);

    expect(exposureSpy).not.toHaveBeenCalled();
    expect(peaksSpy).not.toHaveBeenCalled();
    expect(indicesSpy).not.toHaveBeenCalled();
    expect(groupsSpy).toHaveBeenCalledWith(100);
    expect(groupsSpy).toHaveBeenCalledTimes(1);
  });

  it("multiple exposures: only the cold ones get fetched", async () => {
    const qc = makeQc();
    // Member 100 fully warm.
    qc.setQueryData(queryKeys.exposure(100), { id: 100 });
    qc.setQueryData(queryKeys.peaks(100), []);
    qc.setQueryData(queryKeys.indices(100), []);
    qc.setQueryData(queryKeys.groups(100), []);
    // Member 200 fully cold.
    // Member 300 partially warm (only peaks missing).
    qc.setQueryData(queryKeys.exposure(300), { id: 300 });
    qc.setQueryData(queryKeys.indices(300), []);
    qc.setQueryData(queryKeys.groups(300), []);

    await prefetchColdMembers([100, 200, 300], qc);

    // Exposure: only 200 was missing.
    expect(exposureSpy).toHaveBeenCalledTimes(1);
    expect(exposureSpy).toHaveBeenCalledWith(200);
    // Peaks: 200 + 300 missing.
    expect(peaksSpy).toHaveBeenCalledTimes(2);
    expect(peaksSpy).toHaveBeenCalledWith(200);
    expect(peaksSpy).toHaveBeenCalledWith(300);
    // Indices: only 200 missing.
    expect(indicesSpy).toHaveBeenCalledTimes(1);
    expect(indicesSpy).toHaveBeenCalledWith(200);
    // Groups: only 200 missing.
    expect(groupsSpy).toHaveBeenCalledTimes(1);
    expect(groupsSpy).toHaveBeenCalledWith(200);
  });

  it("awaits all fetches before resolving", async () => {
    // Pin: prefetchColdMembers must not resolve before the underlying
    // fetchQuery promises do — handleSave / buildOverwritePayload need
    // the cache populated before computeMemberSnapshot reads it.
    let resolved = false;
    let resolvePeaks: (v: unknown[]) => void = () => {};
    peaksSpy.mockImplementation(() =>
      new Promise<unknown[]>((res) => { resolvePeaks = res; })
    );
    const qc = makeQc();
    const p = prefetchColdMembers([100], qc).then(() => { resolved = true; });

    // Yield two ticks — the helper has dispatched its inner fetches but
    // they haven't resolved.
    await Promise.resolve(); await Promise.resolve();
    expect(resolved).toBe(false);

    resolvePeaks([]);
    await p;
    expect(resolved).toBe(true);
  });
});
