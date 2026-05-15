import { describe, it, expect } from "vitest";
import { stripQueueMetadata } from "../../src/lib/queue/queueMeta";

describe("stripQueueMetadata", () => {
  it("splits response into meta + payload", () => {
    const response = {
      event_id: 42,
      view_row_id: 99,
      analysis_inputs_hash: "abc",
      client_op_id: "op-1",
      id: 7,
      q: 0.12,
      source: "manual" as const,
    };
    const { meta, payload } = stripQueueMetadata(response);
    expect(meta).toEqual({
      event_id: 42,
      view_row_id: 99,
      analysis_inputs_hash: "abc",
      client_op_id: "op-1",
    });
    expect(payload).toEqual({ id: 7, q: 0.12, source: "manual" });
  });

  it("leaves undefined meta fields as undefined", () => {
    const response = { id: 7 };
    const { meta, payload } = stripQueueMetadata(response);
    expect(meta).toEqual({
      event_id: undefined,
      view_row_id: undefined,
      analysis_inputs_hash: undefined,
      client_op_id: undefined,
    });
    expect(payload).toEqual({ id: 7 });
  });

  it("does not mutate the input", () => {
    const response = { event_id: 1, id: 2 };
    stripQueueMetadata(response);
    expect(response).toEqual({ event_id: 1, id: 2 });
  });
});
