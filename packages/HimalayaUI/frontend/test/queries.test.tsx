import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useExposures, useTrace, usePeaks, useIndices,
  useAddPeak, useRemovePeak, useReanalyzeExposure,
  queryKeys,
} from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  const hasBody = status !== 204 && body !== null;
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(hasBody ? JSON.stringify(body) : null, {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useExposures fetches when sampleId is provided", async () => {
    mockOnce(200, [{ id: 1, sample_id: 10, filename: "f", kind: "file",
                     selected: false, tags: [], sources: [] }]);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useExposures(10), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toHaveLength(1);
  });

  it("useExposures is disabled when sampleId is undefined", () => {
    const fetchSpy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useExposures(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("useTrace fetches for a given exposureId", async () => {
    mockOnce(200, { q: [0.1], I: [10], sigma: [1] });
    const { wrapper } = withClient();
    const { result } = renderHook(() => useTrace(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data?.q).toEqual([0.1]);
  });

  it("usePeaks fetches for a given exposureId", async () => {
    mockOnce(200, []);
    const { wrapper } = withClient();
    const { result } = renderHook(() => usePeaks(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
  });

  it("useIndices fetches for a given exposureId", async () => {
    mockOnce(200, []);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useIndices(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
  });

  it("useAddPeak invalidates peaks and indices for the exposure on success", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.peaks(42), []);
    client.setQueryData(queryKeys.indices(42), []);
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(201, { id: 1, exposure_id: 42, q: 0.1, intensity: null,
                    prominence: null, sharpness: null, source: "manual",
                    stale_indices: 0 });
    const { result } = renderHook(() => useAddPeak(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(0.1); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });

  it("useRemovePeak invalidates peaks and indices for the exposure on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(204, null);
    const { result } = renderHook(() => useRemovePeak(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(7); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });

  it("useReanalyzeExposure invalidates peaks and indices on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 42, analyzed: true });
    const { result } = renderHook(() => useReanalyzeExposure(42), { wrapper });
    await act(async () => { await result.current.mutateAsync(); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(42) });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.indices(42) });
  });
});
