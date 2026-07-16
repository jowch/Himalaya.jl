import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useSeriesTraces } from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(JSON.stringify(body), {
      status,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — useSeriesTraces", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("fetches the exposure_id → Trace map for a series", async () => {
    mockOnce(200, { 1000: { q: [0.1], I: [10], sigma: [1] } });
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSeriesTraces(7), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    // JSON object keys are strings at runtime; number-index access resolves them,
    // exactly as toWaterfallRows does (tracesById[member.exposure_id]).
    expect(result.current.data?.[1000]?.q).toEqual([0.1]);
  });

  it("is disabled (does not fetch) when seriesId is undefined", () => {
    const fetchSpy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSeriesTraces(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});
