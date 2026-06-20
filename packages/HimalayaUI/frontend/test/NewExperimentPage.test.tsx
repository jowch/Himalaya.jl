// test/NewExperimentPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { NewExperimentPage } from "../src/print/pages/NewExperimentPage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><NewExperimentPage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("NewExperimentPage (Phase E1)", () => {
  beforeEach(() => { navigate.mockClear(); });
  afterEach(() => vi.restoreAllMocks());

  it("renders the directory picker", () => {
    renderPage();
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("creates the experiment and routes to its corpus on submit", async () => {
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 682, scanned: 700, message: null });
    vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as api.Experiment);
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/Volumes/data/ssrl/2026_04/1p7m" },
    });
    await waitFor(() => expect(screen.getByTestId("create-submit")).not.toBeDisabled());
    fireEvent.click(screen.getByTestId("create-submit"));
    await waitFor(() => expect(api.createExperiment).toHaveBeenCalled());
    expect(navigate).toHaveBeenCalledWith("/experiments/9/corpus");
  });

  it("keeps submit disabled until validation is ok", () => {
    renderPage();
    expect(screen.getByTestId("create-submit")).toBeDisabled();
  });

  it("blocks submit and warns when the directory is already an experiment", async () => {
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 682, scanned: 700, message: null });
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      { id: 3, name: "Existing", data_dir: "/Volumes/data/ssrl/2026_04/1p7m" } as api.Experiment,
    ]);
    const createSpy = vi.spyOn(api, "createExperiment");
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/Volumes/data/ssrl/2026_04/1p7m" },
    });
    await waitFor(() =>
      expect(screen.getByText(/already (an experiment|uses this directory)/i)).toBeInTheDocument(),
    );
    expect(screen.getByTestId("create-submit")).toBeDisabled();
    fireEvent.click(screen.getByTestId("create-submit"));
    expect(createSpy).not.toHaveBeenCalled();
  });
});
