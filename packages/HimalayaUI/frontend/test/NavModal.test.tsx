import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { renderWithProviders } from "./test-utils";
import { NavModal } from "../src/print/shell/NavModal";
import { useAppState } from "../src/state";
import * as api from "../src/api";

/**
 * NavModal calls useNavigate (M1: commitSample is a door into /sample/:id), so
 * it needs Router context. In production it is a global overlay in CorpusShell,
 * a sibling of the routed content — mirror that here: NavModal beside a tiny
 * Routes whose /sample/:id stub lets a navigation be observed.
 */
function renderModal() {
  return renderWithProviders(
    <MemoryRouter initialEntries={["/samples"]}>
      <NavModal />
      <Routes>
        <Route path="/samples" element={<div data-testid="route-samples" />} />
        <Route
          path="/sample/:sampleId"
          element={<div data-testid="focus-stub" />}
        />
      </Routes>
    </MemoryRouter>,
  );
}

const EXPERIMENTS: api.Experiment[] = [
  { id: 1, name: "SSRL May 2026", path: "/data/ssrl_2026_05", data_dir: "/data/ssrl_2026_05/data", analysis_dir: "/data/ssrl_2026_05/analysis", manifest_path: null, created_at: "2026-05-01", q_units: null, beam_center_x: null, beam_center_y: null, pixel_size_um: null, energy_kev: null, flight_path_m: null },
  { id: 2, name: "APS Apr 2026",  path: "/data/aps_2026_04",  data_dir: "/data/aps_2026_04/data",  analysis_dir: "/data/aps_2026_04/analysis",  manifest_path: null, created_at: "2026-04-15", q_units: null, beam_center_x: null, beam_center_y: null, pixel_size_um: null, energy_kev: null, flight_path_m: null },
];

const SAMPLES_EXP1: api.Sample[] = [
  { id: 10, experiment_id: 1, display_name: "D1", name: "cubic_run03", notes: null, tags: [] },
  { id: 11, experiment_id: 1, display_name: "D2", name: "hex_run01",   notes: null, tags: [] },
];

const SAMPLES_EXP2: api.Sample[] = [
  { id: 20, experiment_id: 2, display_name: "S1", name: "lamellar_A",  notes: null, tags: [] },
];

// Corpus-wide sample list (SA-F4): the experiment step matches sample names
// across ALL experiments through this, not just the committed one.
const CORPUS: api.CorpusSample[] = [
  { id: 10, experiment_id: 1, display_name: "D1",    name: "cubic_run03", notes: null, tags: [], q_units: "A-1" },
  { id: 11, experiment_id: 1, display_name: "D2",    name: "hex_run01",   notes: null, tags: [], q_units: "A-1" },
  { id: 20, experiment_id: 2, display_name: "S1",    name: "lamellar_A",  notes: null, tags: [], q_units: "A-1" },
  { id: 30, experiment_id: 2, display_name: "Calib", name: "aps_calib",   notes: null, tags: [], q_units: "A-1" },
];

function resetStore(): void {
  localStorage.clear();
  useAppState.setState({
    username: "tester",
    activeExperimentId: undefined,
    activeSampleId: undefined,
    activeExposureId: undefined,
    navModalOpen: false,
    navModalStep: "experiment",
  });
}

beforeEach(() => {
  resetStore();
  vi.spyOn(api, "listExperiments").mockResolvedValue(EXPERIMENTS);
  vi.spyOn(api, "listSamples").mockImplementation((expId: number) =>
    Promise.resolve(expId === 1 ? SAMPLES_EXP1 : SAMPLES_EXP2));
  vi.spyOn(api, "listCorpusSamples").mockResolvedValue(CORPUS);
});

describe("<NavModal>", () => {
  it("returns null when navModalOpen is false", () => {
    useAppState.setState({ navModalOpen: false });
    const { container } = renderModal();
    expect(container.querySelector('[data-testid="nav-modal"]')).toBeNull();
  });

  it("renders experiment list when opened at experiment step", async () => {
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderModal();
    expect(await screen.findByText("SSRL May 2026")).toBeInTheDocument();
    expect(screen.getByText("APS Apr 2026")).toBeInTheDocument();
  });

  it("gives each full-bleed result row the inset focus ring so the scroll list can't clip it (UI-RINGCLIP)", async () => {
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderModal();
    await screen.findByText("SSRL May 2026");
    expect(screen.getByTestId("nav-item-experiment-1")).toHaveAttribute(
      "data-focus-ring",
      "inset",
    );
  });

  it("filters experiments by query", async () => {
    const user = userEvent.setup();
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderModal();
    await screen.findByText("SSRL May 2026");

    await user.type(screen.getByTestId("nav-modal-input"), "APS");
    expect(screen.queryByText("SSRL May 2026")).not.toBeInTheDocument();
    // Scoped to the experiment row: the same name can also appear as quiet
    // context on a direct sample hit (SA-F4).
    expect(
      within(screen.getByTestId("nav-item-experiment-2")).getByText("APS Apr 2026"),
    ).toBeInTheDocument();
  });

  it("Enter commits the selected experiment and advances to sample step with chip", async () => {
    const user = userEvent.setup();
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderModal();
    await screen.findByText("SSRL May 2026");

    await user.keyboard("{Enter}"); // commits first item
    await waitFor(() => {
      expect(screen.getByTestId("nav-chip-experiment")).toHaveTextContent("SSRL May 2026");
    });
    await screen.findByText("cubic_run03");
    expect(screen.getByText("hex_run01")).toBeInTheDocument();
  });

  it("Enter on a sample commits the sample and closes the modal", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderModal();
    await screen.findByText("cubic_run03");

    await user.keyboard("{Enter}");
    await waitFor(() => {
      expect(useAppState.getState().navModalOpen).toBe(false);
      expect(useAppState.getState().activeSampleId).toBe(10);
      expect(useAppState.getState().activeExperimentId).toBe(1);
    });
  });

  // M1 on-ramp: committing a sample is the third door into the indexing
  // workspace (beside the contact-sheet status cell and the loupe link). Before
  // this NavModal only wrote the store and left the URL put. cubic_run03 = id 10.
  it("navigates to the focus workspace when a sample is committed", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderModal();
    await screen.findByText("cubic_run03");

    await user.keyboard("{Enter}");
    expect(await screen.findByTestId("focus-stub")).toBeInTheDocument();
  });

  it("Backspace on empty input at sample step pops the experiment chip", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderModal();
    await screen.findByText("cubic_run03");
    expect(screen.getByTestId("nav-chip-experiment")).toBeInTheDocument();

    await user.click(screen.getByTestId("nav-modal-input"));
    await user.keyboard("{Backspace}");

    await waitFor(() => {
      expect(screen.queryByTestId("nav-chip-experiment")).not.toBeInTheDocument();
      expect(screen.getByText("SSRL May 2026")).toBeInTheDocument();
    });
  });

  it("Backspace pops one chip at a time when both are present", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
      activeSampleId: 10,
    });
    renderModal();
    await screen.findByTestId("nav-chip-sample");
    expect(screen.getByTestId("nav-chip-experiment")).toBeInTheDocument();

    await user.click(screen.getByTestId("nav-modal-input"));

    // First backspace — sample chip only
    await user.keyboard("{Backspace}");
    await waitFor(() => {
      expect(screen.queryByTestId("nav-chip-sample")).not.toBeInTheDocument();
    });
    expect(screen.getByTestId("nav-chip-experiment")).toBeInTheDocument();

    // Second backspace — experiment chip
    await user.keyboard("{Backspace}");
    await waitFor(() => {
      expect(screen.queryByTestId("nav-chip-experiment")).not.toBeInTheDocument();
    });
  });

  it("chip × button pops the chip and rewinds", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderModal();
    await screen.findByText("cubic_run03");

    await user.click(screen.getByTestId("nav-chip-experiment-remove"));
    await waitFor(() => {
      expect(screen.queryByTestId("nav-chip-experiment")).not.toBeInTheDocument();
    });
  });

  it("Esc closes without committing", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "experiment",
      activeExperimentId: 99,
      activeSampleId: 42,
    });
    renderModal();
    await screen.findByText("SSRL May 2026");

    await user.click(screen.getByTestId("nav-modal-input"));
    await user.keyboard("{Escape}");
    expect(useAppState.getState().navModalOpen).toBe(false);
    expect(useAppState.getState().activeExperimentId).toBe(99);
    expect(useAppState.getState().activeSampleId).toBe(42);
  });

  // ── SA-F4: direct sample find from the experiment step ──────────────────────
  describe("direct sample find (SA-F4)", () => {
    it("a query matching a sample surfaces the Samples group with experiment context", async () => {
      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      await user.type(screen.getByTestId("nav-modal-input"), "cubic");
      const row = await screen.findByTestId("nav-item-corpus-sample-10");
      expect(screen.getByTestId("nav-samples-group-label")).toHaveTextContent("Samples");
      // Row = sample display name + its experiment name as quiet context.
      expect(within(row).getByText("D1")).toBeInTheDocument();
      expect(within(row).getByText("SSRL May 2026")).toBeInTheDocument();
    });

    it("an empty query shows no Samples group", async () => {
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");
      expect(screen.queryByTestId("nav-samples-group-label")).toBeNull();
    });

    it("an experiments-only query shows no Samples group (no noise)", async () => {
      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      await user.type(screen.getByTestId("nav-modal-input"), "SSRL");
      expect(screen.getByText("SSRL May 2026")).toBeInTheDocument();
      expect(screen.queryByTestId("nav-samples-group-label")).toBeNull();
    });

    it("Enter on a highlighted sample hit navigates to /sample/:id and closes", async () => {
      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      // "cubic" matches no experiment, so the sample hit is the first flat row.
      await user.type(screen.getByTestId("nav-modal-input"), "cubic");
      await screen.findByTestId("nav-item-corpus-sample-10");

      await user.keyboard("{Enter}");
      await waitFor(() => {
        expect(useAppState.getState().navModalOpen).toBe(false);
        expect(useAppState.getState().activeSampleId).toBe(10);
        expect(useAppState.getState().activeExperimentId).toBe(1);
      });
      expect(await screen.findByTestId("focus-stub")).toBeInTheDocument();
    });

    it("clicking a sample hit commits like Enter would", async () => {
      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      await user.type(screen.getByTestId("nav-modal-input"), "hex_run");
      await user.click(await screen.findByTestId("nav-item-corpus-sample-11"));
      await waitFor(() => {
        expect(useAppState.getState().navModalOpen).toBe(false);
        expect(useAppState.getState().activeSampleId).toBe(11);
        expect(useAppState.getState().activeExperimentId).toBe(1);
      });
      expect(await screen.findByTestId("focus-stub")).toBeInTheDocument();
    });

    it("one flat highlight order spans both groups, experiments first", async () => {
      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      // "aps" matches one experiment (APS Apr 2026) and one sample (aps_calib).
      await user.type(screen.getByTestId("nav-modal-input"), "aps");
      const sampleRow = await screen.findByTestId("nav-item-corpus-sample-30");
      const expRow = screen.getByTestId("nav-item-experiment-2");

      // Experiments lead the flat order.
      expect(expRow).toHaveAttribute("data-selected");
      expect(sampleRow).not.toHaveAttribute("data-selected");

      // One ArrowDown crosses the group boundary seamlessly.
      await user.keyboard("{ArrowDown}");
      expect(sampleRow).toHaveAttribute("data-selected");
      expect(expRow).not.toHaveAttribute("data-selected");

      // Enter on the sample hit goes straight to the workspace (exp 2 carried).
      await user.keyboard("{Enter}");
      await waitFor(() => {
        expect(useAppState.getState().navModalOpen).toBe(false);
        expect(useAppState.getState().activeSampleId).toBe(30);
        expect(useAppState.getState().activeExperimentId).toBe(2);
      });
      expect(await screen.findByTestId("focus-stub")).toBeInTheDocument();
    });

    it("caps the group at 8 and discloses the remainder honestly", async () => {
      const many: api.CorpusSample[] = Array.from({ length: 11 }, (_, i) => ({
        id: 100 + i,
        experiment_id: 1,
        display_name: `Bulk ${i + 1}`,
        name: `bulk_${String(i + 1).padStart(2, "0")}`,
        notes: null,
        tags: [],
        q_units: "A-1",
      }));
      vi.spyOn(api, "listCorpusSamples").mockResolvedValue(many);

      const user = userEvent.setup();
      useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
      renderModal();
      await screen.findByText("SSRL May 2026");

      await user.type(screen.getByTestId("nav-modal-input"), "bulk");
      await screen.findByTestId("nav-item-corpus-sample-100");
      const results = screen.getByTestId("nav-modal-results");
      expect(within(results).getAllByTestId(/nav-item-corpus-sample-/)).toHaveLength(8);
      expect(screen.getByTestId("nav-samples-overflow")).toHaveTextContent("+3 more matches");
    });

    it("Backspace rewind semantics are untouched by the Samples group", async () => {
      const user = userEvent.setup();
      useAppState.setState({
        navModalOpen: true,
        navModalStep: "sample",
        activeExperimentId: 1,
      });
      renderModal();
      await screen.findByText("cubic_run03");

      await user.click(screen.getByTestId("nav-modal-input"));
      await user.keyboard("{Backspace}");
      await waitFor(() => {
        expect(screen.queryByTestId("nav-chip-experiment")).not.toBeInTheDocument();
        expect(screen.getByText("SSRL May 2026")).toBeInTheDocument();
      });
      // Rewound to the experiment step with an empty query: no Samples group.
      expect(screen.queryByTestId("nav-samples-group-label")).toBeNull();
    });
  });

  it("clicking a result commits like Enter would", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderModal();
    const row = await screen.findByText("hex_run01");

    await user.click(row);
    await waitFor(() => {
      expect(useAppState.getState().activeSampleId).toBe(11);
      expect(useAppState.getState().navModalOpen).toBe(false);
    });
  });
});
