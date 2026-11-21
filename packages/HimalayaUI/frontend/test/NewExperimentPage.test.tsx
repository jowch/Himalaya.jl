import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { NewExperimentPage } from "../src/print/pages/NewExperimentPage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

const GEOMETRY: api.ManifestGeometry = {
  energy_kev: 9, energy_kev_source: "prp",
  flight_path_m: 1.8095, flight_path_m_source: "setup",
  beam_center_x: 421.4, beam_center_x_source: "setup",
  beam_center_y: 836.9, beam_center_y_source: "setup",
  pixel_size_um: 172, pixel_size_um_source: "computed",
};

function setup(root: string, overrides: Partial<api.SetupPreview> = {}): api.SetupPreview {
  return {
    root,
    name: `20261103/${root.split("/").pop()}`,
    data_dir: `${root}/data`,
    analysis_dir: `${root}/analysis/automatic_analysis/tot_files`,
    setup_file: `${root}/analysis/setup_info_20261103_090000.txt`,
    setup_ambiguous: false,
    layout_detected: true,
    patterns: { image: "{name}_0_001.tif", metadata: "{name}_0_001.prp", integration: "{name}_tot.dat" },
    exposures: 1220,
    without_metadata: 0,
    without_trace: 3,
    geometry: GEOMETRY,
    discrepancies: [],
    duplicate_of: null,
    ...overrides,
  };
}

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><NewExperimentPage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

function pick(path: string): void {
  fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, { target: { value: path } });
}

describe("NewExperimentPage", () => {
  beforeEach(() => {
    navigate.mockClear();
    vi.spyOn(api, "suggestPaths").mockResolvedValue({ suggestions: [] });
  });
  afterEach(() => vi.restoreAllMocks());

  it("puts the cursor in the folder field on arrival", () => {
    renderPage();
    expect(document.activeElement).toBe(screen.getByTestId("dirpicker-input").querySelector("input"));
  });

  it("keeps the add button disabled until a folder with setups is picked", () => {
    renderPage();
    expect(screen.getByTestId("create-submit")).toBeDisabled();
  });

  it("lists every setup under a beamtime folder with a one-line summary", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b/20261103",
      setups: [setup("/b/20261103/setup1"), setup("/b/20261103/setup2", { exposures: 1, without_trace: 0 })],
    });
    renderPage();
    pick("/b/20261103");
    const rows = await screen.findAllByTestId("setup-row");
    expect(rows).toHaveLength(2);
    expect(within(rows[0]!).getByTestId("setup-summary")).toHaveTextContent(
      "1,220 frames · flight path 1.81 m · 3 without a trace yet",
    );
    expect(within(rows[1]!).getByTestId("setup-summary")).toHaveTextContent("1 frame · flight path 1.81 m");
    expect(screen.getByTestId("create-submit")).toHaveTextContent("Add 2 experiments");
  });

  it("creates one experiment per chosen setup with its name, layout and geometry, then goes home", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b/20261103",
      setups: [setup("/b/20261103/setup1"), setup("/b/20261103/setup2")],
    });
    const create = vi.spyOn(api, "createExperiment")
      .mockResolvedValueOnce({ id: 11 } as api.Experiment)
      .mockResolvedValueOnce({ id: 12 } as api.Experiment);
    renderPage();
    pick("/b/20261103");
    const rows = await screen.findAllByTestId("setup-row");
    fireEvent.change(within(rows[1]!).getByTestId("setup-name").querySelector("input")!, { target: { value: "flight path 2" } });
    fireEvent.click(screen.getByTestId("create-submit"));
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments"));
    expect(create).toHaveBeenCalledTimes(2);
    const first = create.mock.calls[0]![0];
    expect(first).toMatchObject({
      name: "20261103/setup1",
      data_dir: "/b/20261103/setup1/data",
      analysis_dir: "/b/20261103/setup1/analysis/automatic_analysis/tot_files",
      patterns: { image: "{name}_0_001.tif", metadata: "{name}_0_001.prp", integration: "{name}_tot.dat" },
    });
    expect(first.geometry).toMatchObject({ flight_path_m: 1.8095, flight_path_m_source: "setup" });
    expect(create.mock.calls[1]![0].name).toBe("flight path 2");
  });

  it("opens the experiment straight away when only one setup is added", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/20261103/setup1", setups: [setup("/b/20261103/setup1")] });
    vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/20261103/setup1");
    await screen.findAllByTestId("setup-row");
    fireEvent.keyDown(window, { key: "Enter", metaKey: true });
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/11"));
  });

  it("⌘/Ctrl+Enter from a focused checkbox or the path field adds the ticked setups and nothing else", async () => {
    vi.mocked(api.suggestPaths).mockResolvedValue({ suggestions: ["/b/20261103/notes"] });
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b/20261103",
      setups: [setup("/b/20261103/setup1"), setup("/b/20261103/setup2")],
    });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/20261103");
    const rows = await screen.findAllByTestId("setup-row");
    const second = within(rows[1]!).getByRole("checkbox");
    fireEvent.click(second);
    fireEvent.keyDown(second, { key: "Enter", ctrlKey: true });
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/11"));
    expect(create).toHaveBeenCalledTimes(1);
    expect(create.mock.calls[0]![0].name).toBe("20261103/setup1");
    expect(second).toHaveAttribute("aria-checked", "false");
  });

  it("⌘/Ctrl+Enter in the path field does not pick a suggestion", async () => {
    vi.mocked(api.suggestPaths).mockResolvedValue({ suggestions: ["/b/20261103/notes"] });
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/20261103", setups: [setup("/b/20261103/setup1")] });
    vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/20261103");
    await screen.findAllByTestId("setup-row");
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input);
    await screen.findByRole("option", { name: "notes" });
    fireEvent.keyDown(input, { key: "Enter", metaKey: true });
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/11"));
    expect(input).toHaveValue("/b/20261103");
  });

  it("plain Enter in the path field adds the found setup rather than entering its first subfolder", async () => {
    vi.mocked(api.suggestPaths).mockResolvedValue({ suggestions: ["/b/20261103/setup1/analysis", "/b/20261103/setup1/data"] });
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/20261103/setup1", setups: [setup("/b/20261103/setup1")] });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/20261103/setup1");
    await screen.findAllByTestId("setup-row");
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input);
    await screen.findByRole("option", { name: "analysis" });
    fireEvent.keyDown(input, { key: "Enter" });
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/11"));
    expect(create).toHaveBeenCalledTimes(1);
    expect(input).toHaveValue("/b/20261103/setup1");
  });

  it("Enter right after picking a subfolder adds only that subfolder's setup, once its own preview has arrived", async () => {
    vi.mocked(api.suggestPaths).mockImplementation(async (p) => ({
      suggestions: p.includes("setup1") ? ["/b/day1/setup1/analysis"] : ["/b/day1/setup1", "/b/day1/setup2"],
    }));
    vi.spyOn(api, "previewFolder").mockImplementation(async (p) =>
      p.includes("setup1")
        ? { path: p, setups: [setup("/b/day1/setup1")] }
        : { path: p, setups: [setup("/b/day1/setup1"), setup("/b/day1/setup2")] },
    );
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/day1/");
    await waitFor(() => expect(screen.getAllByTestId("setup-row")).toHaveLength(2));
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input);
    await screen.findByRole("option", { name: "setup1" });
    fireEvent.keyDown(input, { key: "ArrowDown" });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(input).toHaveValue("/b/day1/setup1/");
    fireEvent.keyDown(input, { key: "Enter" });
    expect(create).not.toHaveBeenCalled();
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/11"));
    expect(create).toHaveBeenCalledTimes(1);
    expect(create.mock.calls[0]![0]).toMatchObject({ data_dir: "/b/day1/setup1/data" });
  });

  it("leaves ⌘/Ctrl+Enter to an open dialog such as the palette", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    const field = document.createElement("input");
    dialog.appendChild(field);
    document.body.appendChild(dialog);
    const reached = vi.fn();
    field.addEventListener("keydown", reached);
    fireEvent.keyDown(field, { key: "Enter", ctrlKey: true });
    expect(reached).toHaveBeenCalled();
    await new Promise((r) => setTimeout(r, 50));
    expect(create).not.toHaveBeenCalled();
    dialog.remove();
  });

  it("adds once when ⌘/Ctrl+Enter is pressed twice quickly", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    let resolve: (e: api.Experiment) => void = () => {};
    const create = vi.spyOn(api, "createExperiment").mockImplementation(() => new Promise((r) => { resolve = r; }));
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    fireEvent.keyDown(window, { key: "Enter", metaKey: true });
    fireEvent.keyDown(window, { key: "Enter", metaKey: true });
    resolve({ id: 11 } as api.Experiment);
    await waitFor(() => expect(navigate).toHaveBeenCalled());
    expect(create).toHaveBeenCalledTimes(1);
  });

  it("leaves a setup whose layout was not recognised unticked", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b", setups: [setup("/b/s1", { layout_detected: false })],
    });
    renderPage();
    pick("/b");
    const rows = await screen.findAllByTestId("setup-row");
    expect(within(rows[0]!).getByRole("checkbox")).toHaveAttribute("aria-checked", "false");
    expect(screen.getByTestId("create-submit")).toBeDisabled();
  });

  it("waits for the recount of an edited setup and sends its geometry", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    let answer: (m: api.ManifestResponse) => void = () => {};
    vi.spyOn(api, "fetchManifest").mockImplementation(() => new Promise((r) => { answer = r; }));
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    fireEvent.click(screen.getByRole("button", { name: "Details" }));
    fireEvent.change(screen.getByTestId("details-setup-file").querySelector("input")!, { target: { value: "/b/s1/other_setup_info.txt" } });
    fireEvent.click(screen.getByTestId("create-submit"));
    await new Promise((r) => setTimeout(r, 20));
    expect(create).not.toHaveBeenCalled();
    answer({ total: 1, matched: { image: 1, metadata: 1, integration: 1 }, unmatched: [],
      geometry: { ...GEOMETRY, flight_path_m: 0.5, flight_path_m_source: "setup" } });
    await waitFor(() => expect(create).toHaveBeenCalled());
    expect(create.mock.calls[0]![0].geometry).toMatchObject({ flight_path_m: 0.5 });
  });

  it("keeps a typed name when the preview refetches", async () => {
    const preview = vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter><NewExperimentPage /></MemoryRouter>
      </QueryClientProvider>,
    );
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    fireEvent.change(screen.getByTestId("setup-name").querySelector("input")!, { target: { value: "kept" } });
    preview.mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1", { exposures: 2000 })] });
    await qc.invalidateQueries();
    await waitFor(() => expect(screen.getByTestId("setup-summary")).toHaveTextContent("2,000 frames"));
    expect(screen.getByTestId("setup-name").querySelector("input")).toHaveValue("kept");
  });

  it("says when a folder cannot be read for another reason", async () => {
    vi.spyOn(api, "previewFolder").mockRejectedValue(new api.ApiError(500, "disk unplugged", null));
    renderPage();
    pick("/b");
    await waitFor(() => expect(screen.getByTestId("preview-status")).toHaveTextContent("Could not read that folder: disk unplugged"));
  });

  it("leaves out a setup that is already an experiment", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b/20261103",
      setups: [setup("/b/20261103/setup1", { duplicate_of: { id: 3, name: "Existing" } }), setup("/b/20261103/setup2")],
    });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 12 } as api.Experiment);
    renderPage();
    pick("/b/20261103");
    const rows = await screen.findAllByTestId("setup-row");
    expect(within(rows[0]!).getByTestId("setup-summary")).toHaveTextContent('Already added as "Existing".');
    expect(screen.getByTestId("create-submit")).toHaveTextContent("Add experiment");
    fireEvent.click(screen.getByTestId("create-submit"));
    await waitFor(() => expect(create).toHaveBeenCalledTimes(1));
    expect(create.mock.calls[0]![0].data_dir).toBe("/b/20261103/setup2/data");
  });

  it("sends a geometry value typed in Details as a user value", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 11 } as api.Experiment);
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    fireEvent.click(screen.getByRole("button", { name: "Details" }));
    fireEvent.change(screen.getByTestId("details-flight_path_m").querySelector("input")!, { target: { value: "2.5" } });
    fireEvent.click(screen.getByTestId("create-submit"));
    await waitFor(() => expect(create).toHaveBeenCalled());
    expect(create.mock.calls[0]![0].geometry).toMatchObject({
      flight_path_m: 2.5, flight_path_m_source: "user", beam_center_x: 421.4, beam_center_x_source: "setup",
    });
  });

  it("recounts from the manifest when a pattern is corrected in Details", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({
      path: "/b/s1",
      setups: [setup("/b/s1", { layout_detected: false, exposures: 0, without_trace: 0, patterns: { image: "{name}.tif", metadata: "{name}.prp", integration: "{name}.dat" } })],
    });
    const manifest = vi.spyOn(api, "fetchManifest").mockResolvedValue({
      total: 10, matched: { image: 5, metadata: 5, integration: 4 },
      unmatched: [{ file: "x", miss: "integration" }], geometry: GEOMETRY,
    });
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    expect(screen.getByTestId("setup-summary")).toHaveTextContent("integration layout not recognised");
    fireEvent.click(screen.getByRole("button", { name: "Details" }));
    fireEvent.change(screen.getByTestId("details-image-pattern").querySelector("input")!, { target: { value: "{name}_0_001.tif" } });
    await waitFor(() => expect(screen.getByTestId("setup-summary")).toHaveTextContent("5 frames"));
    expect(screen.getByTestId("setup-summary")).toHaveTextContent("1 without a trace yet");
    expect(manifest.mock.calls.at(-1)![1]).toMatchObject({ image: "{name}_0_001.tif" });
  });

  it("shows a failed setup's error and stays on the page", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/b/s1", setups: [setup("/b/s1")] });
    vi.spyOn(api, "createExperiment").mockRejectedValue(new Error("an experiment already uses this directory"));
    renderPage();
    pick("/b/s1");
    await screen.findAllByTestId("setup-row");
    fireEvent.click(screen.getByTestId("create-submit"));
    expect(await screen.findByRole("alert")).toHaveTextContent("an experiment already uses this directory");
    expect(navigate).not.toHaveBeenCalled();
  });

  it("says when a folder holds no setups", async () => {
    vi.spyOn(api, "previewFolder").mockResolvedValue({ path: "/empty", setups: [] });
    renderPage();
    pick("/empty");
    await waitFor(() => expect(screen.getByTestId("preview-status")).toHaveTextContent("No images found"));
    expect(screen.getByTestId("create-submit")).toBeDisabled();
  });
});
