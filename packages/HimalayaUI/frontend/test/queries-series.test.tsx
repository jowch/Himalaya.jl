import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useSeriesList, useSeries, queryKeys } from "../src/queries";
import * as api from "../src/api";

function wrapper(client = makeClient()) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
}

describe("series read hooks", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("useSeriesList fetches via api.listSeries under queryKeys.seriesList", async () => {
    const spy = vi.spyOn(api, "listSeries").mockResolvedValue([
      // minimal SeriesSummary
      {
        id: 1, title: "S1", description: null, content_hash: "h",
        created_by: null, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null,
        view_grouping_mode: null, view_show_peak_ticks: null,
        view_show_peak_labels: null, last_event_at: null,
        author_username: null, member_count: 0, member_phases: [],
        member_phase_count: 0, has_stale_members: false,
        ordering_variable: null, spans_experiments: false, experiment_name: null,
      },
    ]);
    const { result } = renderHook(() => useSeriesList(), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(spy).toHaveBeenCalledOnce();
    expect(result.current.data?.[0].title).toBe("S1");
    // key-shape regression: the hook uses the pre-existing seriesList key.
    expect(queryKeys.seriesList).toEqual(["series-list"]);
  });

  it("useSeries(id) fetches via api.getSeries and is gated on a defined id", async () => {
    const spy = vi.spyOn(api, "getSeries").mockResolvedValue({
      id: 7, title: "S7", description: null, content_hash: "",
      created_by: null, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: null, view_show_peak_ticks: null,
      view_show_peak_labels: null, ordering_variable: null,
      order_rule: "ascending", state: "draft", members: [], samples: [],
    });
    const { result } = renderHook(() => useSeries(7), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(spy).toHaveBeenCalledWith(7);
    expect(result.current.data?.id).toBe(7);
  });

  it("useSeries(undefined) is disabled (no fetch)", async () => {
    const spy = vi.spyOn(api, "getSeries").mockResolvedValue({} as never);
    const { result } = renderHook(() => useSeries(undefined), { wrapper: wrapper() });
    expect(result.current.fetchStatus).toBe("idle");
    expect(spy).not.toHaveBeenCalled();
  });
});
