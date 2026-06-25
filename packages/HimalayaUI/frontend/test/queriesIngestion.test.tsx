// test/queriesIngestion.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { queryKeys, useLoads } from "../src/queries";
import * as api from "../src/api";

function wrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return ({ children }: { children: React.ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

describe("ingestion queryKeys + hooks (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("loads key is distinct from samples key, and `undefined`-tolerant like samples", () => {
    expect(queryKeys.loads(7)).toEqual(["experiment", 7, "loads"]);
    expect(queryKeys.samples(7)).toEqual(["experiment", 7, "samples"]);
    // Mirrors queryKeys.samples' `id ?? "none"` shape (queries.ts:51) so a
    // disabled (undefined-id) loads query never prefix-collides with an
    // enabled one.
    expect(queryKeys.loads(undefined)).toEqual(["experiment", "none", "loads"]);
  });

  it("useLoads fetches listLoads", async () => {
    const spy = vi.spyOn(api, "listLoads").mockResolvedValue([]);
    const { result } = renderHook(() => useLoads(7), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(spy).toHaveBeenCalledWith(7);
  });
});
