import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useSaveSeries, useCommitSeriesPlate, useDeleteSeries,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { Series } from "../src/api";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function fullSeries(id: number): Series {
  return {
    id, title: "S", description: null, content_hash: "sha256:x",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: null, order_rule: "manual",
    state: "committed", members: [], samples: [],
  };
}

function mockOnce(status: number, body: unknown): ReturnType<typeof vi.spyOn> {
  return vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — series mutation hooks (I3.5b)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("useSaveSeries with an id PATCHes /api/series/:id and mints a client_op_id", async () => {
    const { wrapper } = withClient();
    const fetchSpy = mockOnce(200, fullSeries(5));
    const { result } = renderHook(() => useSaveSeries(), { wrapper });
    act(() => { result.current.mutate({ id: 5, title: "S", samples: [] }); });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    const init = fetchSpy.mock.calls[0]?.[1];
    expect(url).toBe("/api/series/5");
    expect(init?.method).toBe("PATCH");
    // client_op_id is minted inside mutationFn and rides as a header.
    const headers = new Headers(init?.headers);
    expect(headers.get("X-Client-Op-Id")).toBeTruthy();
  });

  it("useCommitSeriesPlate POSTs /api/series/:id/commit (spinner)", async () => {
    const { wrapper } = withClient();
    const fetchSpy = mockOnce(200, fullSeries(5));
    const { result } = renderHook(() => useCommitSeriesPlate(), { wrapper });
    act(() => { result.current.mutate({ id: 5, members: [] }); });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    expect(url).toBe("/api/series/5/commit");
    expect(fetchSpy.mock.calls[0]?.[1]?.method).toBe("POST");
  });

  it("useDeleteSeries DELETEs /api/series/:id", async () => {
    const { wrapper } = withClient();
    const fetchSpy = mockOnce(200, { id: 5, deleted: true, event_id: 9 });
    const { result } = renderHook(() => useDeleteSeries(), { wrapper });
    act(() => { result.current.mutate({ id: 5 }); });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    expect(url).toBe("/api/series/5");
    expect(fetchSpy.mock.calls[0]?.[1]?.method).toBe("DELETE");
  });
});
