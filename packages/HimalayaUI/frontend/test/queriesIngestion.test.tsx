import { describe, it, expect, vi, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { queryKeys, useProblems } from "../src/queries";
import * as api from "../src/api";

function wrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return ({ children }: { children: React.ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

describe("ingest problems query", () => {
  afterEach(() => vi.restoreAllMocks());

  it("keys problems under the experiment", () => {
    expect(queryKeys.problems(7)).toEqual(["experiment", 7, "problems"]);
  });

  it("useProblems fetches listProblems for the experiment", async () => {
    const spy = vi.spyOn(api, "listProblems").mockResolvedValue({ ingest_error: null, exposures: [] });
    const { result } = renderHook(() => useProblems(7), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(spy).toHaveBeenCalledWith(7);
  });
});
