// test/App.ingestListener.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { PrintApp } from "../src/print/App";
import { useAppState } from "../src/state";
import * as api from "../src/api";

// Minimal EventSource stub that captures listeners so the test can fire frames.
class FakeES {
  static last: FakeES | null = null;
  listeners: Record<string, ((e: MessageEvent) => void)[]> = {};
  url: string;
  constructor(url: string) { this.url = url; FakeES.last = this; }
  addEventListener(t: string, fn: (e: MessageEvent) => void) {
    (this.listeners[t] ??= []).push(fn);
  }
  removeEventListener() {}
  close() {}
  emit(type: string, data: unknown) {
    for (const fn of this.listeners[type] ?? []) fn({ data: JSON.stringify(data) } as MessageEvent);
  }
}

describe("App ingestInFlight SSE listener (Phase E1)", () => {
  beforeEach(() => {
    (globalThis as unknown as { EventSource: unknown }).EventSource = FakeES;
    useAppState.setState({ ingestInFlight: null });
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => { useAppState.setState({ ingestInFlight: null }); vi.restoreAllMocks(); });

  function mount() {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments"]}><PrintApp /></MemoryRouter>
      </QueryClientProvider>,
    );
  }

  it("writes ingestInFlight on ingest_progress", () => {
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_progress", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 300, total: 680 },
    });
    expect(useAppState.getState().ingestInFlight?.[7]).toEqual({
      processed: 300, total: 680, status: "scanning",
    });
  });

  it("maps a rescan-phase frame to status=analyzing (the rescanning surface)", () => {
    // An initial scan (create route) sends no phase → "scanning" → GroupingReviewPage.
    // A rescan (/{id}/scan) tags its frames phase:"rescan" → "analyzing" → the
    // inline ProgressBar, since the experiment's table data is already present.
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_progress", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 120, total: 400, phase: "rescan" },
    });
    expect(useAppState.getState().ingestInFlight?.[7]).toEqual({
      processed: 120, total: 400, status: "analyzing",
    });
  });

  it("a rescan-phase ingest_started also yields status=analyzing", () => {
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_started", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 0, total: 0, phase: "rescan" },
    });
    expect(useAppState.getState().ingestInFlight?.[7]?.status).toBe("analyzing");
  });

  it("clears ingestInFlight on ingest_complete", () => {
    mount();
    FakeES.last!.emit("curation", { kind: "ingest_started", entity_id: 7, entity_type: "experiment", payload: { experiment_id: 7, processed: 0, total: 680 } });
    FakeES.last!.emit("curation", { kind: "ingest_complete", entity_id: 7, entity_type: "experiment", payload: { experiment_id: 7, processed: 680, total: 680 } });
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("an ingest_* frame is short-circuited (does NOT reach handleRemoteEvent)", async () => {
    // Spy on the queue reconciler: an ingest frame must NOT trigger it (it is
    // broadcast-only, never an own-op; running it would roll back pending edits).
    const rc = await import("../src/lib/queue/replayCoordinator");
    const spy = vi.spyOn(rc, "handleRemoteEvent");
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_progress", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 300, total: 680 },
    });
    expect(spy).not.toHaveBeenCalled();
    // A NON-ingest frame still flows through (control).
    FakeES.last!.emit("curation", { kind: "peak_added", entity_id: 1, entity_type: "exposure", payload: {} });
    expect(spy).toHaveBeenCalledTimes(1);
  });
});
