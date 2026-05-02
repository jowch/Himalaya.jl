import { describe, expect, it, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { handleCurationEvent } from "../src/lib/sseSubscriber";

describe("handleCurationEvent", () => {
  it("invalidates peaks/indices/groups/exposure queries for the event's entity_id", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", entity_id: 7, actor: "bob" }),
      { username: "alice", qc },
    );
    expect(invalidate).toHaveBeenCalledTimes(4);
    const keys = invalidate.mock.calls.map((c) => (c[0] as { queryKey: unknown }).queryKey);
    expect(keys).toContainEqual(["exposure", 7, "peaks"]);
    expect(keys).toContainEqual(["exposure", 7, "indices"]);
    expect(keys).toContainEqual(["exposure", 7, "groups"]);
    expect(keys).toContainEqual(["exposure-entity", 7]);
  });

  it("skips self-echo events whose actor matches local username", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", entity_id: 7, actor: "alice" }),
      { username: "alice", qc },
    );
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("ignores malformed JSON without crashing", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent("{not json", { username: "alice", qc });
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("ignores non-exposure entity_types", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "x", entity_type: "experiment", entity_id: 7, actor: null }),
      { username: "alice", qc },
    );
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("propagates anonymous events (actor: null) from other tabs", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", entity_id: 3, actor: null }),
      { username: "alice", qc },
    );
    expect(invalidate).toHaveBeenCalledTimes(4);
  });

  it("ignores events with missing entity_id", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", actor: "bob" }),
      { username: "alice", qc },
    );
    expect(invalidate).not.toHaveBeenCalled();
  });
});
