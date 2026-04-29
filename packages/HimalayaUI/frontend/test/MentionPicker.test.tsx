import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionPicker } from "../src/components/MentionPicker";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const INDICES: api.IndexEntry[] = [
  { id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91, r_squared: 0.998,
    lattice_d: 5.14, status: "candidate", peaks: [], predicted_q: [1.223, 1.414] },
];
const PEAKS: api.Peak[] = [
  { id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
    sharpness: 0.3, source: "auto", excluded: false },
];
const EXPOSURES: api.Exposure[] = [
  { id: 1, sample_id: 3, filename: "JC001-120.dat", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "", tags: [], sources: [] },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeSampleId: 3, activeExposureId: 1, activeExperimentId: 1 });
  vi.spyOn(api, "listIndices").mockResolvedValue(INDICES);
  vi.spyOn(api, "listPeaks").mockResolvedValue(PEAKS);
  vi.spyOn(api, "listExposures").mockResolvedValue(EXPOSURES);
  vi.spyOn(api, "listSamples").mockResolvedValue([]);
});

describe("<MentionPicker>", () => {
  it("shows indices section when query matches phase", async () => {
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/Pn3m/)).toBeInTheDocument();
  });

  it("shows peak section when query matches q value", async () => {
    renderWithProviders(
      <MentionPicker query="1.22" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/1\.223/)).toBeInTheDocument();
  });

  it("calls onSelect with [[type:id]] token when row is clicked", async () => {
    const onSelect = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={onSelect} onDismiss={vi.fn()} />
    );
    await user.click(await screen.findByText(/Pn3m/));
    expect(onSelect).toHaveBeenCalledWith("[[index:17]]");
  });

  it("calls onDismiss on Escape key", async () => {
    const onDismiss = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={onDismiss} />
    );
    await user.keyboard("{Escape}");
    expect(onDismiss).toHaveBeenCalled();
  });

  it("shows empty state when no results match", async () => {
    renderWithProviders(
      <MentionPicker query="zzznothing" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/no results/i)).toBeInTheDocument();
  });
});
