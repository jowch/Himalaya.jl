import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useDeleteIndex,
  queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { IndexEntry } from "../src/api";

const EXPOSURE_ID = 42;

function makeIndex(id: number): IndexEntry {
  return {
    id,
    exposure_id: EXPOSURE_ID,
    phase: "Pn3m",
    basis: 100,
    score: 0.9,
    r_squared: 0.99,
    lattice_d: 50,
    ngc: null,
    status: "candidate",
    kind: "speculative",
    inputs_hash: null,
    peaks: [],
    predicted_q: [],
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

describe("queries — index deletion (queue-driven, M2.3)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  // ---------------------- useDeleteIndex ----------------------

  it("useDeleteIndex filters the index from the indices cache optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    mockNever();
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
      expect(indices.map((i) => i.id)).toEqual([102]);
    });
  });

  it("useDeleteIndex success leaves the optimistic state in place", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    mockOnce(200, { deleted: 101 });
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
    expect(indices.map((i) => i.id)).toEqual([102]);
  });

  it("useDeleteIndex rolls back the indices cache on error", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    mockOnce(403, { error: "only speculative indices can be deleted" });
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
      expect(indices.map((i) => i.id)).toEqual([101, 102]);
    });
  });
});
