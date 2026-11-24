import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, act } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter, useLocation } from "react-router";
import { App } from "../src/App";
import * as api from "../src/api";
import * as toast from "../src/lib/toast";
import { getClientId } from "../src/lib/clientId";
import { FakeEventSource } from "./helpers/fakeEventSource";

function Where() {
  return <output aria-label="where">{useLocation().pathname}</output>;
}

function mount(path: string): QueryClient {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><App /><Where /></MemoryRouter>
    </QueryClientProvider>,
  );
  return qc;
}

function deleted(id: number, clientId: string | null = "another-window") {
  act(() => {
    FakeEventSource.current().emit("curation", {
      kind: "experiment_deleted", entity_type: "experiment", entity_id: id,
      payload: null, client_id: clientId, client_op_id: null, id: 1,
    });
  });
}

describe("an experiment removed in another window", () => {
  beforeEach(() => {
    FakeEventSource.install();
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockRejectedValue(new Error("404"));
  });
  afterEach(() => { vi.restoreAllMocks(); });

  it("drops its cached data and refreshes the experiment list and folder previews", () => {
    const qc = mount("/experiments");
    qc.setQueryData(["experiment", 7], { id: 7 });
    qc.setQueryData(["experiment", 7, "samples"], []);
    qc.setQueryData(["experiment", 8], { id: 8 });
    const spy = vi.spyOn(qc, "invalidateQueries");
    deleted(7);
    expect(qc.getQueryData(["experiment", 7])).toBeUndefined();
    expect(qc.getQueryData(["experiment", 7, "samples"])).toBeUndefined();
    expect(qc.getQueryData(["experiment", 8])).toEqual({ id: 8 });
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(["experiments"]));
    expect(keys).toContain(JSON.stringify(["folder-preview"]));
  });

  it("sends a window that has it open back home and says why", () => {
    const toastSpy = vi.spyOn(toast, "showToast");
    mount("/experiments/7/workspace/3");
    deleted(7);
    expect(screen.getByLabelText("where").textContent).toMatch(/^\/experiments$/);
    expect(toastSpy).toHaveBeenCalledWith("This experiment was removed in another window.");
  });

  it("leaves a window on another experiment where it is", () => {
    const toastSpy = vi.spyOn(toast, "showToast");
    mount("/experiments/70/workspace/3");
    deleted(7);
    expect(screen.getByLabelText("where").textContent).toBe("/experiments/70/workspace/3");
    expect(toastSpy).not.toHaveBeenCalled();
  });

  it("does not move or notify the window that removed it", () => {
    const toastSpy = vi.spyOn(toast, "showToast");
    mount("/experiments/7/workspace/3");
    deleted(7, getClientId());
    expect(screen.getByLabelText("where").textContent).toBe("/experiments/7/workspace/3");
    expect(toastSpy).not.toHaveBeenCalled();
  });
});
