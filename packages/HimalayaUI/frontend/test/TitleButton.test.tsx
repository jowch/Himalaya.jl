import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { TitleButton } from "../src/components/TitleButton";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const EXPERIMENT: api.Experiment = {
  id: 1, name: "SSRL May 2026", path: "/data/x", data_dir: "/data/x/data",
  analysis_dir: "/data/x/analysis", manifest_path: null, created_at: "2026-05-01",
};

const SAMPLES: api.Sample[] = [
  { id: 10, experiment_id: 1, label: "D1", name: "cubic_run03", notes: null, tags: [] },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    username: "tester",
    activeExperimentId: undefined,
    activeSampleId: undefined,
    navModalOpen: false,
    navModalStep: "experiment",
  });
  vi.spyOn(api, "getExperiment").mockResolvedValue(EXPERIMENT);
  vi.spyOn(api, "listSamples").mockResolvedValue(SAMPLES);
});

describe("<TitleButton>", () => {
  it("renders the empty-state prompt when nothing is selected", () => {
    renderWithProviders(<TitleButton />);
    expect(screen.getByText(/pick an experiment/i)).toBeInTheDocument();
    expect(screen.getByTestId("title-button-kbd")).toHaveTextContent("/");
  });

  it("renders the sample name as the title and the experiment as detail", async () => {
    useAppState.setState({ activeExperimentId: 1, activeSampleId: 10 });
    renderWithProviders(<TitleButton />);
    expect(await screen.findByText("cubic_run03")).toBeInTheDocument();
    expect(await screen.findByText("SSRL May 2026")).toBeInTheDocument();
  });

  it("opens the nav modal at the experiment step when no experiment set", async () => {
    const user = userEvent.setup();
    renderWithProviders(<TitleButton />);
    await user.click(screen.getByTestId("title-button"));
    expect(useAppState.getState().navModalOpen).toBe(true);
    expect(useAppState.getState().navModalStep).toBe("experiment");
  });

  it("opens the nav modal at the sample step when experiment is already set", async () => {
    const user = userEvent.setup();
    useAppState.setState({ activeExperimentId: 1 });
    renderWithProviders(<TitleButton />);
    await user.click(screen.getByTestId("title-button"));
    expect(useAppState.getState().navModalOpen).toBe(true);
    expect(useAppState.getState().navModalStep).toBe("sample");
  });
});
