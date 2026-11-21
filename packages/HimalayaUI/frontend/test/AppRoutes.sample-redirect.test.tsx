import { describe, it, expect, vi, afterEach } from "vitest";
import { StrictMode } from "react";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as toast from "../src/lib/toast";

vi.mock("../src/workspace/WorkspacePage", () => ({ WorkspacePage: () => <div data-testid="workspace-stub" /> }));

function Where(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="where">{loc.pathname + loc.search}</div>;
}

function stubApi(resolveBody?: unknown): void {
  vi.spyOn(global, "fetch").mockImplementation((input) => {
    const url = String(input);
    const body = url.includes("/api/resolve")
      ? resolveBody
      : url.endsWith("/api/samples")
        ? [{ id: 42, experiment_id: 7, name: "A1", notes: null, q_units: "A-1" }]
        : [];
    return Promise.resolve({ ok: true, status: 200, json: () => Promise.resolve(body), text: () => Promise.resolve(JSON.stringify(body ?? null)) } as Response);
  });
}

function renderAt(route: string): void {
  render(
    <StrictMode>
      <QueryClientProvider client={new QueryClient({ defaultOptions: { queries: { retry: false } } })}>
        <MemoryRouter initialEntries={[route]}>
          <AppRoutes />
          <Where />
        </MemoryRouter>
      </QueryClientProvider>
    </StrictMode>,
  );
}

function serverDown(): void {
  vi.spyOn(global, "fetch").mockImplementation(() => Promise.reject(new TypeError("Failed to fetch")));
}

afterEach(() => vi.restoreAllMocks());

describe("old sample links open the workspace", () => {
  it.each(["/sample/42", "/sample/42/loupe"])("%s lands on the sample in its experiment's workspace", async (route) => {
    stubApi();
    renderAt(route);
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments/7/workspace/42"));
  });

  it("a sample that no longer exists lands on the experiments home and says the link is stale", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    stubApi();
    renderAt("/sample/99");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments"));
    expect(shown).toHaveBeenCalledTimes(1);
    expect(shown).toHaveBeenCalledWith("That sample link is no longer valid", "warning");
  });

  it("an old sample link opened while the server is down says it couldn't open, not that it is stale", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    serverDown();
    renderAt("/sample/42");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments"));
    expect(shown).toHaveBeenCalledTimes(1);
    expect(shown).toHaveBeenCalledWith("Couldn't open that sample link", "error");
  });

  it("a permalink opened while the server is down says it couldn't open, not that it is stale", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    serverDown();
    renderAt("/index/lipid/A1");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments"));
    expect(shown).toHaveBeenCalledTimes(1);
    expect(shown).toHaveBeenCalledWith("Couldn't open that sample link", "error");
  });

  it("a permalink that no longer resolves lands on the experiments home and says the link is stale", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    stubApi({ error: "not found" });
    renderAt("/index/lipid/gone");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments"));
    expect(shown).toHaveBeenCalledWith("That sample link is no longer valid", "warning");
  });

  it("a link that resolves says nothing", async () => {
    const shown = vi.spyOn(toast, "showToast").mockImplementation(() => {});
    stubApi({ experiment_id: 7, experiment_name: "lipid", sample_id: 42, sample_name: "A1" });
    renderAt("/index/lipid/A1");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments/7/workspace/42"));
    expect(shown).not.toHaveBeenCalled();
  });

  it("an old /index/:experiment/:sample permalink resolves the slugs and opens the workspace", async () => {
    stubApi({ experiment_id: 7, experiment_name: "lipid", sample_id: 42, sample_name: "A1" });
    renderAt("/index/lipid/A1");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments/7/workspace/42"));
  });

  it("keeps the frame an old Loupe link asked for", async () => {
    stubApi();
    renderAt("/sample/42/loupe?experiment=7&exposure=105");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments/7/workspace/42?exposure=105"));
  });

  it("keeps the frame a resolved permalink names", async () => {
    stubApi({ experiment_id: 7, experiment_name: "lipid", sample_id: 42, sample_name: "A1", exposure_id: 105 });
    renderAt("/index/lipid/A1");
    await waitFor(() => expect(screen.getByTestId("where").textContent).toBe("/experiments/7/workspace/42?exposure=105"));
  });
});
