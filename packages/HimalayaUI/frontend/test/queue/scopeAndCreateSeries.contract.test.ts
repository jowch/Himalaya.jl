import { describe, it, expect, vi, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { scopeAndCreateSeriesMutator } from "../../src/lib/queue/mutators/scopeAndCreateSeries";
import * as api from "../../src/api";

// Full `Series` fixture — matches the real GET/POST projection shape
// (api.isFullSeries true: id + state + members[] + samples[]).
const FULL_SERIES: api.Series = {
  id: 7,
  title: "Series by ratio",
  description: null,
  content_hash: "abc",
  created_by: null,
  created_at: null,
  updated_at: null,
  forked_from_id: null,
  forked_at_hash: null,
  forked_from_title: null,
  view_grouping_mode: null,
  view_show_peak_ticks: null,
  view_show_peak_labels: null,
  ordering_variable: "ratio",
  order_rule: "ascending",
  state: "draft",
  members: [],
  samples: [],
};

describe("scopeAndCreateSeriesMutator", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("kind is series_save (reuses the existing series SSE event kind)", () => {
    expect(scopeAndCreateSeriesMutator.kind).toBe("series_save");
  });

  it("onMutate is a no-op (user navigates away on success; nothing to patch)", () => {
    const qc = new QueryClient();
    const ctx = scopeAndCreateSeriesMutator.onMutate(
      {
        key: "ratio",
        tags: [],
        title: "S",
        samples: [],
        username: "a",
        clientId: "c",
        clientOpId: "op1",
      } as never,
      qc,
    );
    expect(typeof ctx.restore).toBe("function");
    expect(() => ctx.restore()).not.toThrow();
  });

  it("request: calls batchSampleTags (source='scoping') THEN saveSeries (create, no id)", async () => {
    const order: string[] = [];
    const batchSpy = vi
      .spyOn(api, "batchSampleTags")
      .mockImplementation(async () => {
        order.push("batch");
        return [];
      });
    const saveSpy = vi.spyOn(api, "saveSeries").mockImplementation(async () => {
      order.push("save");
      return FULL_SERIES;
    });
    const result = await scopeAndCreateSeriesMutator.request(
      {
        key: "ratio",
        tags: [{ sampleId: 10, value: "1:1" }],
        title: "Series by ratio",
        samples: [{ sample_id: 10, position: 0 }],
        orderingVariable: "ratio",
        username: "a",
        clientId: "c",
        clientOpId: "op1",
      } as never,
      new AbortController().signal,
    );
    // Tags are written before the series is created.
    expect(order).toEqual(["batch", "save"]);
    expect(batchSpy).toHaveBeenCalledWith(
      "ratio",
      [{ sample_id: 10, value: "1:1" }],
      "scoping",
      expect.anything(),
    );
    expect(saveSpy).toHaveBeenCalledWith(
      expect.objectContaining({
        title: "Series by ratio",
        ordering_variable: "ratio",
        samples: [{ sample_id: 10, position: 0 }],
      }),
      undefined, // no id → POST /api/series (create)
      expect.anything(),
    );
    // request returns the created Series.
    expect(result).toBe(FULL_SERIES);
  });

  it("onSuccess: splices the new series into the series-list cache and invalidates scoping caches", () => {
    const qc = new QueryClient();
    qc.setQueryData(["series-list"], []);
    const inv = vi.spyOn(qc, "invalidateQueries");
    scopeAndCreateSeriesMutator.onSuccess(
      {
        key: "ratio",
        tags: [],
        title: "S",
        samples: [],
        username: "a",
        clientId: "c",
        clientOpId: "op1",
      } as never,
      FULL_SERIES,
      qc,
    );
    // Splices the full series into the series-list.
    expect(qc.getQueryData(["series-list"])).toEqual([FULL_SERIES]);
    // Primes the detail cache so the /series/:id navigation reads it immediately.
    expect(qc.getQueryData(["series", 7])).toEqual(FULL_SERIES);
    // Invalidates the scoping proposal sources.
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
  });

  it("defines no synthesizeFromSse (series_save foreign frames route to saveSeriesMutator)", () => {
    expect(scopeAndCreateSeriesMutator.synthesizeFromSse).toBeUndefined();
  });
});
