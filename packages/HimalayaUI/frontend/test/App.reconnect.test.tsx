import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { PrintApp } from "../src/print/App";
import * as api from "../src/api";
import { FakeEventSource } from "./helpers/fakeEventSource";

const KEY = ["samples", "corpus-probe"];
const OTHER = ["exposure", 9, "peaks"];

describe("PrintApp after the live-update stream reconnects", () => {
  let qc: QueryClient;

  beforeEach(() => {
    FakeEventSource.install();
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    qc = new QueryClient({ defaultOptions: { queries: { retry: false, staleTime: Infinity } } });
    qc.setQueryData(KEY, ["cached before the backend restarted"]);
    qc.setQueryData(OTHER, []);
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments"]}><PrintApp /></MemoryRouter>
      </QueryClientProvider>,
    );
  });
  afterEach(() => vi.restoreAllMocks());

  it("keeps the cache on the first connection", () => {
    FakeEventSource.current().fire("open");
    expect(qc.getQueryState(KEY)?.isInvalidated).toBe(false);
  });

  it("marks every cached query stale when the stream opens again, since frames sent while it was down are lost", () => {
    FakeEventSource.current().fire("open");
    FakeEventSource.current().fire("error");
    FakeEventSource.current().fire("open");
    expect(qc.getQueryState(KEY)?.isInvalidated).toBe(true);
    expect(qc.getQueryState(OTHER)?.isInvalidated).toBe(true);
  });
});
