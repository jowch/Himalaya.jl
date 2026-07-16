/**
 * `treats404AsSuccess` contract.
 *
 * Idempotent removes (peakRemove, deleteIndex,
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
import {
  useRemovePeak, useDeleteIndex, useRemoveSampleTag,
  queryKeys,
} from "../../src/queries";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import { deleteIndexMutator } from "../../src/lib/queue/mutators/indexGroup";
import {
  removeSampleTagMutator, removeExposureTagMutator,
} from "../../src/lib/queue/mutators/trivial";
import { setToastImpl } from "../../src/lib/toast";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import type { Peak, IndexEntry, Sample } from "../../src/api";

describe("treats404AsSuccess flag — set on idempotent remove mutators", () => {
  it("peakRemove", () => {
    expect(peakRemoveMutator.treats404AsSuccess).toBe(true);
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

// The framework branch is mutator-shape-agnostic in principle, but a regression
// dependent on rollback-context shape, response type, or cache topology would
// only be caught by exercising more than one mutator. Two extra cases — across
// distinct cache shapes (indices, samples-with-nested-tags) — pin this.
describe("treats404AsSuccess framework branch — coverage across mutator shapes", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    setToastImpl(() => {});
  });
  afterEach(() => { setToastImpl(null); });

  function mockFetch404(): void {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "missing" }), {
        status: 404, headers: { "Content-Type": "application/json" },
      }),
    );
  }

  it("useDeleteIndex: optimistic removal sticks in the indices cache on 404", async () => {
    const EXPOSURE_ID = 5;
    const INDEX: IndexEntry = {
      id: 10, exposure_id: EXPOSURE_ID, phase: "Pn3m", basis: 0.1, score: 0.9,
      r_squared: 0.99, lattice_d: 50, ngc: 0.5, status: "candidate",
      kind: "auto", inputs_hash: "h1", peaks: [], predicted_q: [0.1],
    };
    const client = makeClient();
    client.setQueryData(queryKeys.indices(EXPOSURE_ID), [INDEX]);
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    mockFetch404();
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(10); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    expect(client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID))).toEqual([]);
  });

  it("useRemoveSampleTag: optimistic tag removal sticks on 404", async () => {
    const EXPERIMENT_ID = 1;
    const SAMPLE_ID = 10;
    const SAMPLE: Sample = {
      id: SAMPLE_ID, experiment_id: EXPERIMENT_ID, name: "n",
      notes: null,
      tags: [{ id: 99, key: "buffer", value: "PBS", source: "manual" }],
    };
    const client = makeClient();
    client.setQueryData(queryKeys.samples(EXPERIMENT_ID), [SAMPLE]);
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    mockFetch404();
    const { result } = renderHook(
      () => useRemoveSampleTag(EXPERIMENT_ID, SAMPLE_ID), { wrapper },
    );
    act(() => { result.current.mutate(99); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    const sample = client.getQueryData<Sample[]>(queryKeys.samples(EXPERIMENT_ID))?.[0];
    expect(sample?.tags).toEqual([]);
  });
});
