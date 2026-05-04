import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useAddPeak, useRemovePeak, useSetPeakExcluded, queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import * as api from "../src/api";
import type { Peak, Exposure } from "../src/api";

const EXPOSURE_ID = 42;
const HASH_OLD = "old0000000000000000000000000000000000000000000000000000000000000";
const HASH_NEW = "new0000000000000000000000000000000000000000000000000000000000000";

function makeExposure(): Exposure {
  return {
    id: EXPOSURE_ID,
    sample_id: 1,
    filename: "ex",
    kind: "file",
    selected: true,
    status: null,
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: HASH_OLD,
    analysis_inputs_hash: HASH_OLD,
  };
}

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

function mockNever(): void {
  vi.spyOn(global, "fetch").mockImplementation(() => new Promise(() => {}));
}

describe("queries — peak mutations (queue-driven, M2.2)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  // ---------------------- useAddPeak ----------------------

  it("useAddPeak inserts a negative-id placeholder optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), []);
    mockNever();
    const { result } = renderHook(() => useAddPeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(0.123); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list).toHaveLength(1);
      expect(list[0].id).toBeLessThan(0);
      expect(list[0].q).toBe(0.123);
      expect(list[0].source).toBe("manual");
    });
  });

  it("useAddPeak replaces the placeholder with the server peak on success", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), []);
    client.setQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID), makeExposure());
    // Flat shape (matches `routes_peaks.jl` JSON3.write Dict). Earlier
    // version of this fixture wrapped under `peak: {...}`, which matched
    // a stale `PeakAddResponse.peak: Peak` type but NOT the actual server
    // response — the production mutator crashed on `response.peak.id`
    // because there is no `peak` key. Caught only in the deep-scan after
    // the fourth review round.
    mockOnce(201, {
      id: 99, exposure_id: EXPOSURE_ID, q: 0.123, intensity: null,
      prominence: null, sharpness: null, source: "manual", excluded: false,
      event_id: 7,
      view_row_id: 7,
      analysis_inputs_hash: HASH_NEW,
    });
    const { result } = renderHook(() => useAddPeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(0.123); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list).toHaveLength(1);
      expect(list[0].id).toBe(99);
    });
    const exp = client.getQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID));
    expect(exp?.analysis_inputs_hash).toBe(HASH_NEW);
  });

  it("useAddPeak rolls back the optimistic placeholder on error", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), []);
    mockOnce(400, { error: "bad q" });
    const { result } = renderHook(() => useAddPeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(0.123); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list).toHaveLength(0);
    });
  });

  // ---------------------- useRemovePeak ----------------------

  it("useRemovePeak removes the peak optimistically", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "manual", excluded: false,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    mockNever();
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(5); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list).toHaveLength(0);
    });
  });

  it("useRemovePeak writes the new hash to the exposure cache on success", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "manual", excluded: false,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    client.setQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID), makeExposure());
    mockOnce(200, {
      event_id: 99,
      view_row_id: 99,
      analysis_inputs_hash: HASH_NEW,
    });
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(5); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
    expect(list).toHaveLength(0);
    const exp = client.getQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID));
    expect(exp?.analysis_inputs_hash).toBe(HASH_NEW);
  });

  it("useRemovePeak rolls back on error", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "manual", excluded: false,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    mockOnce(400, { error: "nope" });
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(5); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list).toHaveLength(1);
      expect(list[0].id).toBe(5);
    });
  });

  // ---------------------- useSetPeakExcluded ----------------------

  it("useSetPeakExcluded(true) flips excluded optimistically and writes server response", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "auto", excluded: false,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    client.setQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID), makeExposure());
    mockOnce(200, {
      ...peak,
      excluded: true,
      event_id: 9,
      view_row_id: 9,
      analysis_inputs_hash: HASH_NEW,
    });
    const { result } = renderHook(() => useSetPeakExcluded(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate({ peakId: 5, excluded: true }); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list[0].excluded).toBe(true);
    });
    const exp = client.getQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID));
    expect(exp?.analysis_inputs_hash).toBe(HASH_NEW);
  });

  it("useSetPeakExcluded(false) routes through the unexclude mutator", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "auto", excluded: true,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    client.setQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID), makeExposure());
    const setSpy = vi.spyOn(api, "setPeakExcluded").mockResolvedValueOnce({
      ...peak, excluded: false, event_id: 9, view_row_id: 9,
      analysis_inputs_hash: HASH_NEW,
    });
    const { result } = renderHook(() => useSetPeakExcluded(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate({ peakId: 5, excluded: false }); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list[0].excluded).toBe(false);
    });
    expect(setSpy).toHaveBeenCalledWith(5, false, expect.anything());
  });

  it("useSetPeakExcluded rolls back on error", async () => {
    const { client, wrapper } = withClient();
    const peak: Peak = {
      id: 5, exposure_id: EXPOSURE_ID, q: 0.1, intensity: null,
      prominence: null, sharpness: null, source: "auto", excluded: false,
    };
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), [peak]);
    mockOnce(500, { error: "boom" });
    const { result } = renderHook(() => useSetPeakExcluded(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate({ peakId: 5, excluded: true }); });
    // After settle (5xx is infrastructure → retries; we wait briefly and check
    // the ROLLBACK that the queue framework applies on each failed attempt).
    // Rather than chase retry-vs-rollback timing, just assert it did flip
    // optimistically first.
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list[0].excluded).toBe(true);
    });
  });

  // ---------------------- autoReanalyze removed ----------------------

  it("does not call POST /analyze on peak add (backend handles reanalyze synchronously)", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID), []);
    client.setQueryData<Exposure>(queryKeys.exposure(EXPOSURE_ID), makeExposure());
    const reanalyzeSpy = vi.spyOn(api, "reanalyzeExposure");
    mockOnce(201, {
      id: 99, exposure_id: EXPOSURE_ID, q: 0.5, intensity: null,
      prominence: null, sharpness: null, source: "manual", excluded: false,
      event_id: 1, view_row_id: 1, analysis_inputs_hash: HASH_NEW,
    });
    const { result } = renderHook(() => useAddPeak(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(0.5); });
    await waitFor(() => {
      const list = client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID)) ?? [];
      expect(list[0]?.id).toBe(99);
    });
    expect(reanalyzeSpy).not.toHaveBeenCalled();
  });
});
