import { describe, expect, it, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { handleCurationEvent } from "../src/lib/sseSubscriber";
import { queryKeys } from "../src/queries";

describe("handleCurationEvent", () => {
  it("invalidates peaks/indices/groups/exposure queries for the event's entity_id", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", entity_id: 7, actor: "bob" }),
      { clientId: "test-tab", qc },
    );
    expect(invalidate).toHaveBeenCalledTimes(4);
    const keys = invalidate.mock.calls.map((c) => (c[0] as { queryKey: unknown }).queryKey);
    expect(keys).toContainEqual(["exposure", 7, "peaks"]);
    expect(keys).toContainEqual(["exposure", 7, "indices"]);
    expect(keys).toContainEqual(["exposure", 7, "groups"]);
    expect(keys).toContainEqual(["exposure-entity", 7]);
  });

  it("ignores malformed JSON without crashing", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent("{not json", { clientId: "test-tab", qc });
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("ignores non-exposure entity_types", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "x", entity_type: "experiment", entity_id: 7, actor: null }),
      { clientId: "test-tab", qc },
    );
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("propagates anonymous events (actor: null) from other tabs", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", entity_id: 3, actor: null }),
      { clientId: "test-tab", qc },
    );
    expect(invalidate).toHaveBeenCalledTimes(4);
  });

  it("ignores events with missing entity_id", () => {
    const invalidate = vi.fn();
    const qc = { invalidateQueries: invalidate } as unknown as QueryClient;
    handleCurationEvent(
      JSON.stringify({ kind: "peak_added", entity_type: "exposure", actor: "bob" }),
      { clientId: "test-tab", qc },
    );
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("skips events authored by this tab (client_id match)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    handleCurationEvent(
      JSON.stringify({
        kind: "peak_added",
        entity_type: "exposure",
        entity_id: 5,
        actor: "alice",
        client_id: "tab-A",
      }),
      { clientId: "tab-A", qc },
    );
    expect(spy).not.toHaveBeenCalled();
  });

  it("processes events from another tab even when same actor", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    handleCurationEvent(
      JSON.stringify({
        kind: "peak_added",
        entity_type: "exposure",
        entity_id: 5,
        actor: "alice",
        client_id: "tab-B",
      }),
      { clientId: "tab-A", qc },
    );
    expect(spy).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(5) });
  });

  it("processes system events (NULL client_id) regardless of tab", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    handleCurationEvent(
      JSON.stringify({
        kind: "analyze_run",
        entity_type: "exposure",
        entity_id: 5,
        actor: null,
        client_id: null,
      }),
      { clientId: "tab-A", qc },
    );
    expect(spy).toHaveBeenCalled();
  });

  it("ignores malformed JSON", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    handleCurationEvent("{not json", { clientId: "tab-A", qc });
    expect(spy).not.toHaveBeenCalled();
  });
});
