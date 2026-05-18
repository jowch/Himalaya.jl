/**
 * applyRemoteToCache — series foreign-event handlers (#166/#167/#168).
 *
 * series_created / series_recipe_updated carry no post_state (master plan
 * §5.2) and the SSE payload's recipe rows are id-less + replay-volatile (§11),
 * so the handler is invalidate-only. series_deleted removes the detail cache
 * and filters the listing.
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";

describe("applyRemoteToCache: series #166 kinds", () => {
  it("series_created invalidates the detail + listing caches", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    qc.setQueryData(queryKeys.seriesList, [{ id: 5 }]);
    const remote: SseEvent = {
      id: 1, kind: "series_created", entity_type: "series", entity_id: 5,
      payload: { title: "S", samples: [] },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryState(queryKeys.series(5))?.isInvalidated).toBe(true);
    expect(qc.getQueryState(queryKeys.seriesList)?.isInvalidated).toBe(true);
  });

  it("series_recipe_updated invalidates the detail cache", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    const remote: SseEvent = {
      id: 2, kind: "series_recipe_updated", entity_type: "series", entity_id: 5,
      payload: { samples: [] },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryState(queryKeys.series(5))?.isInvalidated).toBe(true);
  });

  it("series_deleted removes the detail cache and filters the listing", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(5), { id: 5 });
    qc.setQueryData(queryKeys.seriesList, [{ id: 5 }, { id: 6 }]);
    const remote: SseEvent = {
      id: 3, kind: "series_deleted", entity_type: "series", entity_id: 5,
      payload: { id: 5 },
    };
    applyRemoteToCache(remote, qc);
    expect(qc.getQueryData(queryKeys.series(5))).toBeUndefined();
    expect(qc.getQueryData(queryKeys.seriesList)).toEqual([{ id: 6 }]);
  });
});
