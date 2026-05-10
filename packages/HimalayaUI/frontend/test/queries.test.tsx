import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useExposures, useTrace, usePeaks, useIndices,
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

  // Regression for #116. The hook used to take an `excludeRejected` option
  // that landed in the cache key, splitting the same sample's exposures
  // across two cache rows (Index page filtered, Inspect/MentionPicker not).
  // Result: cold re-fetch on every Index↔Inspect crossing. Pinning the key
  // shape keeps a future regressor from quietly reintroducing the suffix.
  it("useExposures keys on the bare exposures prefix (no opts in cache key)", async () => {
    mockOnce(200, []);
    const { client, wrapper } = withClient();
    const { result } = renderHook(() => useExposures(10), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(client.getQueryData(["sample", 10, "exposures"] as const))
      .toEqual([]);
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

  // useAddPeak / useRemovePeak / useSetPeakExcluded migrated to the queue
  // framework (M2.2) — see test/queries-peaks.test.tsx for the new coverage.
  // useReanalyzeExposure migrated to the queue framework (M2.5) — see
  // test/queries-reanalyze.test.tsx for the new coverage.
});
