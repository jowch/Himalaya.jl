/**
 * peakRemove invariant: optimistic placeholder ids (negative) must never reach
 * the mutator. The TraceViewer click handler is responsible for skipping
 * negative-id peaks; the mutator's onMutate throws as defense-in-depth so any
 * regression at the click layer fails loudly instead of producing the silent
 * "ghost peak" desync (DELETE /api/peaks/-N → 404 → rollback re-inserts the
 * placeholder → original add resolves into a real, undeletable-feeling peak).
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import { queryKeys } from "../../src/queries";
import type { Peak } from "../../src/api";

function makeQc(): QueryClient {
  return new QueryClient({ defaultOptions: { queries: { retry: false } } });
}

const PEAK: Peak = {
  id: 7, exposure_id: 5, q: 0.5, intensity: 1.0, prominence: 0.8,
  sharpness: 30, source: "manual", excluded: false,
};

describe("peakRemoveMutator invariant", () => {
  it("throws when called with a negative (optimistic placeholder) id", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.peaks(5), [PEAK]);
    expect(() => {
      peakRemoveMutator.onMutate({
        kind: "peak_removed", clientOpId: "op",
        exposureId: 5, username: "alice", clientId: "tab",
        peakId: -3,
      } as never, qc);
    }).toThrow(/optimistic placeholder/i);
  });

  it("does not mutate the cache when the invariant fires", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.peaks(5), [PEAK]);
    const before = JSON.parse(JSON.stringify(qc.getQueryData(queryKeys.peaks(5))));
    try {
      peakRemoveMutator.onMutate({
        kind: "peak_removed", clientOpId: "op",
        exposureId: 5, username: "alice", clientId: "tab",
        peakId: -3,
      } as never, qc);
    } catch { /* expected */ }
    expect(qc.getQueryData(queryKeys.peaks(5))).toEqual(before);
  });

  it("works normally for positive (server-confirmed) ids", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.peaks(5), [PEAK]);
    const ctx = peakRemoveMutator.onMutate({
      kind: "peak_removed", clientOpId: "op",
      exposureId: 5, username: "alice", clientId: "tab",
      peakId: 7,
    } as never, qc);
    expect(qc.getQueryData<Peak[]>(queryKeys.peaks(5))).toEqual([]);
    ctx.restore();
    expect(qc.getQueryData<Peak[]>(queryKeys.peaks(5))).toEqual([PEAK]);
  });
});
