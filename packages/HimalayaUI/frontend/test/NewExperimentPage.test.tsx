// test/NewExperimentPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { NewExperimentPage } from "../src/print/pages/NewExperimentPage";
import { useDraftExperiment } from "../src/lib/draftExperiment";
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
  beforeEach(() => {
    navigate.mockClear();
    useDraftExperiment.getState().clear();
  });
  afterEach(() => vi.restoreAllMocks());

  it("renders the directory picker", () => {
    renderPage();
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("Review commits the path to a draft and navigates to draft Configuration WITHOUT creating", async () => {
    const create = vi.spyOn(api, "createExperiment");
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 5, scanned: 5, message: null });
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/data/run42" },
    });
    await waitFor(() => expect(screen.getByRole("button", { name: /review/i })).not.toBeDisabled());
    fireEvent.click(screen.getByRole("button", { name: /review/i }));
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/new/config"));
    expect(create).not.toHaveBeenCalled();
    // The picker commits the experiment ROOT; Configuration resolves it.
    expect(useDraftExperiment.getState().root).toBe("/data/run42");
  });

  it("keeps submit disabled until validation is ok", () => {
    renderPage();
    expect(screen.getByRole("button", { name: /review/i })).toBeDisabled();
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
    expect(screen.getByRole("button", { name: /review/i })).toBeDisabled();
    fireEvent.click(screen.getByRole("button", { name: /review/i }));
    expect(createSpy).not.toHaveBeenCalled();
  });

  it("shows the structure check when the resolver finds a data/analysis layout", async () => {
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 5, scanned: 5, message: null });
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "mini-1p7m", data_dir: "/d/data", analysis_dir: "/d/analysis",
      setup_file: "/d/analysis/setup.txt", setup_ambiguous: false,
      image_pattern: "{name}_0_001.tif", metadata_pattern: "{name}_0_001.prp",
      integration_pattern: "{name}_tot.dat",
    });
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/d" },
    });
    await waitFor(() =>
      expect(screen.getByText(/looks like an experiment/i)).toBeInTheDocument(),
    );
  });

  it("structure check is advisory (does not block) when no layout is detected", async () => {
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 1, scanned: 1, message: null });
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "raw", data_dir: "/raw", analysis_dir: null,
      setup_file: null, setup_ambiguous: true,
      image_pattern: null, metadata_pattern: null, integration_pattern: null,
    });
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/raw" },
    });
    // The advisory note appears as a non-error status …
    const note = await screen.findByText(/no data\/analysis layout detected/i);
    expect(note).toHaveAttribute("role", "status");
    // … and Review stays enabled (exists + not-dup still pass).
    await waitFor(() => expect(screen.getByRole("button", { name: /review/i })).not.toBeDisabled());
  });

  it("Cancel clears the draft and navigates to /experiments", async () => {
    renderPage();
    fireEvent.click(screen.getByRole("button", { name: /cancel/i }));
    expect(navigate).toHaveBeenCalledWith("/experiments");
    expect(useDraftExperiment.getState().root).toBe("");
  });
});
