import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { PrintApp } from "../src/print/App";
import * as api from "../src/api";

class FakeES {
  static last: FakeES | null = null;
  listeners: Record<string, ((e: Event) => void)[]> = {};
  constructor(public url: string) { FakeES.last = this; }
  addEventListener(t: string, fn: (e: Event) => void) {
    (this.listeners[t] ??= []).push(fn);
  }
  removeEventListener() {}
  close() {}
  fire(type: string) {
    for (const fn of this.listeners[type] ?? []) fn(new Event(type));
  }
}

const KEY = ["samples", "corpus-probe"];

describe("PrintApp after the live-update stream reconnects", () => {
  let qc: QueryClient;

  beforeEach(() => {
    (globalThis as unknown as { EventSource: unknown }).EventSource = FakeES;
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    qc = new QueryClient({ defaultOptions: { queries: { retry: false, staleTime: Infinity } } });
    qc.setQueryData(KEY, ["cached before the backend restarted"]);
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments"]}><PrintApp /></MemoryRouter>
      </QueryClientProvider>,
    );
  });
  afterEach(() => vi.restoreAllMocks());

  it("keeps the cache on the first connection", () => {
    FakeES.last!.fire("open");
    expect(qc.getQueryState(KEY)?.isInvalidated).toBe(false);
  });

  it("marks every cached query stale when the stream opens again, since frames sent while it was down are lost", () => {
    FakeES.last!.fire("open");
    FakeES.last!.fire("error");
    FakeES.last!.fire("open");
    expect(qc.getQueryState(KEY)?.isInvalidated).toBe(true);
  });
});
