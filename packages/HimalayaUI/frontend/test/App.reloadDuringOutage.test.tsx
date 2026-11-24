import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { StrictMode } from "react";
import { render, act, cleanup } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router";
import { App } from "../src/App";
import * as api from "../src/api";
import * as toast from "../src/lib/toast";
import { STORAGE_KEY, SCHEMA_VERSION } from "../src/lib/queue/persistence";
import { FakeEventSource } from "./helpers/fakeEventSource";

describe("App mounting after a reload during an outage", () => {
  let qc: QueryClient;

  beforeEach(() => {
    FakeEventSource.install();
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    sessionStorage.clear();
  });
  afterEach(() => {
    cleanup();
    qc.clear();
    vi.restoreAllMocks();
  });

  it("replays the queued edit once, under StrictMode, and keeps it in storage until it lands", async () => {
    const put = vi.spyOn(api, "setExposureStatus").mockRejectedValue(
      new api.NetworkError("PUT", "/api/exposures/7/status"),
    );
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([{
      schemaVersion: SCHEMA_VERSION,
      kind: "set_exposure_status",
      clientOpId: "op-before-reload",
      payload: {
        kind: "set_exposure_status",
        clientOpId: "op-before-reload",
        exposureId: 7,
        sampleId: 3,
        status: "rejected",
      },
    }]));
    qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    await act(async () => {
      render(
        <StrictMode>
          <QueryClientProvider client={qc}>
            <MemoryRouter initialEntries={["/experiments"]}><App /></MemoryRouter>
          </QueryClientProvider>
        </StrictMode>,
      );
    });
    expect(put).toHaveBeenCalled();
    const stored = JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]");
    expect(stored.map((o: { clientOpId: string }) => o.clientOpId)).toEqual(["op-before-reload"]);
    expect(qc.getMutationCache().getAll().filter((m) => m.state.status === "pending")).toHaveLength(1);
  });

  it("tells the user when a saved edit can't be replayed", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([null]));
    qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    await act(async () => {
      render(
        <QueryClientProvider client={qc}>
          <MemoryRouter initialEntries={["/experiments"]}><App /></MemoryRouter>
        </QueryClientProvider>,
      );
    });
    expect(shown).toHaveBeenCalledWith("1 edits from a previous session couldn't be restored", "warning");
  });
});
