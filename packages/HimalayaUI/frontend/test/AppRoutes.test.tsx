import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { MemoryRouter, useNavigate } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes, PAGE_GONE } from "../src/shell/AppRoutes";
import * as api from "../src/api";
import * as toast from "../src/lib/toast";

function GoTo({ path }: { path: string }): JSX.Element {
  const navigate = useNavigate();
  return <button type="button" onClick={() => navigate(path)}>go</button>;
}

function renderRoutes(initialPath: string, goTo?: string) {
  const qc = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: Infinity, staleTime: 0 }, mutations: { retry: false } },
  });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <AppRoutes />
        {goTo && <GoTo path={goTo} />}
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes", () => {
  beforeEach(() => {
    useAppState.setState({ paletteOpen: false });
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("bare / lands on the experiments home inside the app shell", async () => {
    renderRoutes("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
    expect(screen.getByTestId("app-shell")).toBeInTheDocument();
  });

  it.each(["/totally/unknown", "/series/4", "/sample/3", "/index/lipid/JC404", "/samples", "/inspect", "/experiments/7/corpus", "/experiments/7/grouping", "/experiments/new/config"])(
    "sends %s home with a toast",
    async (path) => {
      const shown = vi.spyOn(toast, "showToast");
      renderRoutes(path);
      expect(await screen.findByText("All experiments")).toBeInTheDocument();
      expect(shown).toHaveBeenCalledWith(PAGE_GONE, "warning");
    },
  );

  it("mounts the workspace outside the app shell", async () => {
    vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockResolvedValue({ id: 7, name: "BL-19" } as Awaited<ReturnType<typeof api.getExperiment>>);
    renderRoutes("/experiments/7/workspace");
    expect(await screen.findByTestId("workspace")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
    expect(await screen.findByText("No samples in this experiment.")).toBeInTheDocument();
  });

  it("starts a fresh workspace when the experiment changes, so nothing carries over", async () => {
    vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockImplementation(async (id: number) => ({ id, name: `E${id}` }) as Awaited<ReturnType<typeof api.getExperiment>>);
    renderRoutes("/experiments/7/workspace", "/experiments/8/workspace");
    const first = await screen.findByTestId("workspace");
    fireEvent.click(screen.getByRole("button", { name: "go" }));
    const second = await screen.findByTestId("workspace");
    expect(second).not.toBe(first);
  });
});
