import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { saveSeriesMutator } from "../../src/lib/queue/mutators/saveSeries";
import { deleteSeriesMutator } from "../../src/lib/queue/mutators/deleteSeries";
import { commitSeriesPlateMutator } from "../../src/lib/queue/mutators/commitSeriesPlate";
import { queryKeys } from "../../src/queries";
import type { Series } from "../../src/api";

function fullSeries(id: number): Series {
  return {
    id, title: "S", description: null, content_hash: "sha256:abc",
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, ordering_variable: null, order_rule: "manual",
    state: "draft", members: [], samples: [],
  };
}

describe("saveSeriesMutator", () => {
  it("onSuccess writes a full Series response into the detail cache", () => {
    const qc = new QueryClient();
    const resp = fullSeries(7);
    // `as any` for the flat-payload arg — the established test convention
    // (see deleteComparison.test.tsx); a precise FlatPayload literal is noise.
    saveSeriesMutator.onSuccess(
      { kind: "series_save", payload: {}, clientOpId: "op1",
        title: "S", samples: [], username: "alice", clientId: "c1" } as any,
      resp, qc);
    expect(qc.getQueryData<Series>(queryKeys.series(7))).toEqual(resp);
  });

  it("synthesizeFromSse yields a partial shape carrying the entity id", () => {
    const synth = saveSeriesMutator.synthesizeFromSse?.(
      { id: 99, kind: "series_created", entity_type: "series", entity_id: 7,
        payload: { title: "S", samples: [] } },
      { event_id: 99, client_op_id: "op1", analysis_inputs_hash: undefined });
    expect(synth?.id).toBe(7);
  });
});

describe("deleteSeriesMutator", () => {
  it("onSuccess removes the detail cache and filters the listing", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.series(7), fullSeries(7));
    qc.setQueryData(queryKeys.seriesList, [{ id: 7 }, { id: 8 }]);
    deleteSeriesMutator.onSuccess(
      { kind: "series_delete", payload: {}, clientOpId: "op1",
        id: 7, username: "alice", clientId: "c1" } as any,
      { id: 7, deleted: true, event_id: 99 }, qc);
    expect(qc.getQueryData(queryKeys.series(7))).toBeUndefined();
    expect(qc.getQueryData(queryKeys.seriesList)).toEqual([{ id: 8 }]);
  });
});

describe("commitSeriesPlateMutator", () => {
  it("synthesizeFromSse builds a full Series from the post_state envelope", () => {
    const post = fullSeries(7);
    post.state = "committed";
    const synth = commitSeriesPlateMutator.synthesizeFromSse?.(
      { id: 99, kind: "series_plate_committed", entity_type: "series",
        entity_id: 7, payload: { members: [] }, post_state: post },
      { event_id: 99, client_op_id: "op1", analysis_inputs_hash: undefined });
    expect(synth?.state).toBe("committed");
    expect(synth?.id).toBe(7);
  });
});
