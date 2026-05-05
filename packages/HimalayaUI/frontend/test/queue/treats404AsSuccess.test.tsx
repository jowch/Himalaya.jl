/**
 * `treats404AsSuccess` contract.
 *
 * Idempotent removes (peakRemove, removeIndexFromGroup, deleteIndex,
 * removeSampleTag, removeExposureTag) treat HTTP 404 as a no-op success:
 * the optimistic effect already reflects "row gone", and a 404 means the
 * server has it gone too — restoring would re-insert a phantom row visible
 * until the next refetch.
 *
 * Failure mode being prevented: a 5xx-then-retry where the first attempt
 * succeeded (server deleted the row + emitted SSE) but the response was lost
 * to a network blip. The retry hits a 404, classified as a validation error,
 * which without this flag triggers rollback → cache shows a row the server
 * (and SSE) have already removed.
 */
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "../test-utils";
import { useRemovePeak } from "../../src/queries";
import { queryKeys } from "../../src/queries";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  removeIndexFromGroupMutator, deleteIndexMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import {
  removeSampleTagMutator, removeExposureTagMutator,
} from "../../src/lib/queue/mutators/trivial";
import { setToastImpl } from "../../src/lib/toast";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import type { Peak } from "../../src/api";

describe("treats404AsSuccess flag — set on idempotent remove mutators", () => {
  it("peakRemove", () => {
    expect(peakRemoveMutator.treats404AsSuccess).toBe(true);
  });
  it("removeIndexFromGroup", () => {
    expect(removeIndexFromGroupMutator.treats404AsSuccess).toBe(true);
  });
  it("deleteIndex", () => {
    expect(deleteIndexMutator.treats404AsSuccess).toBe(true);
  });
  it("removeSampleTag", () => {
    expect(removeSampleTagMutator.treats404AsSuccess).toBe(true);
  });
  it("removeExposureTag", () => {
    expect(removeExposureTagMutator.treats404AsSuccess).toBe(true);
  });
});

describe("treats404AsSuccess framework branch — useRemovePeak under HTTP 404", () => {
  const EXPOSURE_ID = 5;
  const PEAK: Peak = {
    id: 7, exposure_id: EXPOSURE_ID, q: 0.5, intensity: 1.0,
    prominence: 0.8, sharpness: 30, source: "manual", excluded: false,
  };

  let toastCalls: Array<{ msg: string; kind: string }> = [];

  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    toastCalls = [];
    setToastImpl((msg, kind) => { toastCalls.push({ msg, kind }); });
  });
  afterEach(() => {
    setToastImpl(null);
  });

  function withSeededClient() {
    const client = makeClient();
    client.setQueryData(queryKeys.peaks(EXPOSURE_ID), [PEAK]);
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    return { client, wrapper };
  }

  function mockFetch404(): void {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "peak not found" }), {
        status: 404, headers: { "Content-Type": "application/json" },
      }),
    );
  }

  it("keeps the optimistic delete in cache (does not rollback)", async () => {
    const { client, wrapper } = withSeededClient();
    mockFetch404();
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(PEAK.id); });
    await waitFor(() =>
      expect(client.getQueryState(queryKeys.peaks(EXPOSURE_ID))?.fetchStatus).toBe("idle"),
    );
    // The peak should be GONE from the cache — optimistic filter held, no
    // rollback fired. Without the flag, restore() would re-insert PEAK.
    expect(client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID))).toEqual([]);
  });

  it("does not surface a validation toast", async () => {
    const { wrapper } = withSeededClient();
    mockFetch404();
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(PEAK.id); });
    // Wait long enough for onError to have run if it were going to.
    await waitFor(() => expect(result.current.isPending).toBe(false));
    expect(toastCalls).toEqual([]);
  });
});
