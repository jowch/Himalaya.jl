import { describe, it, expect, vi, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { createSeriesMutator } from "../../src/lib/queue/mutators/createSeries";
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

describe("createSeriesMutator", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("kind is series_save (the create route emits a series_created wire frame)", () => {
    expect(createSeriesMutator.kind).toBe("series_save");
  });

  it("onMutate is a no-op (user navigates away on success; nothing to patch)", () => {
    const qc = new QueryClient();
    const ctx = createSeriesMutator.onMutate(
      { title: "S", samples: [], username: "a", clientId: "c", clientOpId: "op1" } as never,
      qc,
    );
    expect(typeof ctx.restore).toBe("function");
    expect(() => ctx.restore()).not.toThrow();
  });

  it("request: a SINGLE write — saveSeries with no id (POST create)", async () => {
    const saveSpy = vi.spyOn(api, "saveSeries").mockResolvedValue(FULL_SERIES);
    const batchSpy = vi.spyOn(api, "batchSampleTags").mockResolvedValue([]);
    const result = await createSeriesMutator.request(
      {
        title: "Series by ratio",
        samples: [{ sample_id: 10, position: 0 }],
        ordering_variable: "ratio",
        username: "a",
        clientId: "c",
        clientOpId: "op1",
      } as never,
      new AbortController().signal,
    );
    // ONLY the create runs — no tag write inside this mutator (that is Op A).
    expect(batchSpy).not.toHaveBeenCalled();
    expect(saveSpy).toHaveBeenCalledTimes(1);
    expect(saveSpy).toHaveBeenCalledWith(
      expect.objectContaining({
        title: "Series by ratio",
        ordering_variable: "ratio",
        samples: [{ sample_id: 10, position: 0 }],
      }),
      undefined, // no id → POST /api/series (create)
      expect.objectContaining({ clientOpId: "op1" }),
    );
    // request returns the created Series — reliable because it's a single write.
    expect(result).toBe(FULL_SERIES);
  });

  it("onSuccess: splices the new series into the series-list cache, primes detail, invalidates scoping caches", () => {
    const qc = new QueryClient();
    qc.setQueryData(["series-list"], []);
    const inv = vi.spyOn(qc, "invalidateQueries");
    createSeriesMutator.onSuccess(
      { title: "S", samples: [], username: "a", clientId: "c", clientOpId: "op1" } as never,
      FULL_SERIES,
      qc,
    );
    // Splices the full series into the series-list.
    expect(qc.getQueryData(["series-list"])).toEqual([FULL_SERIES]);
    // Primes the detail cache so the /series/:id navigation reads it immediately.
    expect(qc.getQueryData(["series", 7])).toEqual(FULL_SERIES);
    // Invalidates the scoping proposal sources (the preceding tag write changed them).
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
  });

  it("defines no synthesizeFromSse (series_created foreign frames route to saveSeriesMutator)", () => {
    expect(createSeriesMutator.synthesizeFromSse).toBeUndefined();
  });
});
