import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { NavModal } from "../src/components/NavModal";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const EXPERIMENTS: api.Experiment[] = [
  { id: 1, name: "SSRL May 2026", path: "/data/ssrl_2026_05", data_dir: "/data/ssrl_2026_05/data", analysis_dir: "/data/ssrl_2026_05/analysis", manifest_path: null, created_at: "2026-05-01" },
  { id: 2, name: "APS Apr 2026",  path: "/data/aps_2026_04",  data_dir: "/data/aps_2026_04/data",  analysis_dir: "/data/aps_2026_04/analysis",  manifest_path: null, created_at: "2026-04-15" },
];

const SAMPLES_EXP1: api.Sample[] = [
  { id: 10, experiment_id: 1, label: "D1", name: "cubic_run03", notes: null, tags: [] },
  { id: 11, experiment_id: 1, label: "D2", name: "hex_run01",   notes: null, tags: [] },
];

const SAMPLES_EXP2: api.Sample[] = [
  { id: 20, experiment_id: 2, label: "S1", name: "lamellar_A",  notes: null, tags: [] },
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
});

describe("<NavModal>", () => {
  it("returns null when navModalOpen is false", () => {
    useAppState.setState({ navModalOpen: false });
    const { container } = renderWithProviders(<NavModal />);
    expect(container.querySelector('[data-testid="nav-modal"]')).toBeNull();
  });

  it("renders experiment list when opened at experiment step", async () => {
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderWithProviders(<NavModal />);
    expect(await screen.findByText("SSRL May 2026")).toBeInTheDocument();
    expect(screen.getByText("APS Apr 2026")).toBeInTheDocument();
  });

  it("filters experiments by query", async () => {
    const user = userEvent.setup();
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderWithProviders(<NavModal />);
    await screen.findByText("SSRL May 2026");

    await user.type(screen.getByTestId("nav-modal-input"), "APS");
    expect(screen.queryByText("SSRL May 2026")).not.toBeInTheDocument();
    expect(screen.getByText("APS Apr 2026")).toBeInTheDocument();
  });

  it("Enter commits the selected experiment and advances to sample step with chip", async () => {
    const user = userEvent.setup();
    useAppState.setState({ navModalOpen: true, navModalStep: "experiment" });
    renderWithProviders(<NavModal />);
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
    renderWithProviders(<NavModal />);
    await screen.findByText("cubic_run03");

    await user.keyboard("{Enter}");
    await waitFor(() => {
      expect(useAppState.getState().navModalOpen).toBe(false);
      expect(useAppState.getState().activeSampleId).toBe(10);
      expect(useAppState.getState().activeExperimentId).toBe(1);
    });
  });

  it("Backspace on empty input at sample step pops the experiment chip", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderWithProviders(<NavModal />);
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
    renderWithProviders(<NavModal />);
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
    renderWithProviders(<NavModal />);
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
    renderWithProviders(<NavModal />);
    await screen.findByText("SSRL May 2026");

    await user.click(screen.getByTestId("nav-modal-input"));
    await user.keyboard("{Escape}");
    expect(useAppState.getState().navModalOpen).toBe(false);
    expect(useAppState.getState().activeExperimentId).toBe(99);
    expect(useAppState.getState().activeSampleId).toBe(42);
  });

  it("clicking a result commits like Enter would", async () => {
    const user = userEvent.setup();
    useAppState.setState({
      navModalOpen: true,
      navModalStep: "sample",
      activeExperimentId: 1,
    });
    renderWithProviders(<NavModal />);
    const row = await screen.findByText("hex_run01");

    await user.click(row);
    await waitFor(() => {
      expect(useAppState.getState().activeSampleId).toBe(11);
      expect(useAppState.getState().navModalOpen).toBe(false);
    });
  });
});
