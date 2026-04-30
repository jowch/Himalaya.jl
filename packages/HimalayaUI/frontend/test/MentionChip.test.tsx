import { describe, it, expect, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionChip } from "../src/components/MentionChip";
import { useAppState } from "../src/state";
import type { ResolvedMention } from "../src/hooks/useMentionResolution";
import type { Peak, IndexEntry } from "../src/api";

const PEAK: Peak = {
  id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
  sharpness: 0.3, source: "auto", excluded: false,
};

const INDEX: IndexEntry = {
  id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91,
  r_squared: 0.998, lattice_d: 5.14, ngc: -1.51, status: "candidate",
  peaks: [], predicted_q: [1.223, 1.414, 1.732],
};

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ hoveredPeakId: undefined, hoveredIndexId: undefined });
});

describe("<MentionChip> — resolved", () => {
  it("renders peak chip with q value", () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    renderWithProviders(<MentionChip resolved={resolved} />);
    expect(screen.getByText(/q = 1\.223/)).toBeInTheDocument();
  });

  it("renders index chip with phase and score", () => {
    const resolved: ResolvedMention = { type: "index", data: INDEX };
    renderWithProviders(<MentionChip resolved={resolved} />);
    expect(screen.getByText(/Pn3m/)).toBeInTheDocument();
    expect(screen.getByText(/0\.91/)).toBeInTheDocument();
  });

  it("sets hoveredPeakId on mouse enter for peak chip", async () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/q = 1\.223/));
    expect(useAppState.getState().hoveredPeakId).toBe(42);
  });

  it("clears hoveredPeakId on mouse leave for peak chip", async () => {
    const resolved: ResolvedMention = { type: "peak", data: PEAK };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/q = 1\.223/));
    await user.unhover(screen.getByText(/q = 1\.223/));
    expect(useAppState.getState().hoveredPeakId).toBeUndefined();
  });

  it("sets hoveredIndexId on mouse enter for index chip", async () => {
    const resolved: ResolvedMention = { type: "index", data: INDEX };
    const user = userEvent.setup();
    renderWithProviders(<MentionChip resolved={resolved} />);
    await user.hover(screen.getByText(/Pn3m/));
    expect(useAppState.getState().hoveredIndexId).toBe(17);
  });
});

describe("<MentionChip> — loading and dead", () => {
  it("renders original text for 'loading' state", () => {
    renderWithProviders(<MentionChip resolved="loading" originalText="Pn3m · 0.91" />);
    expect(screen.getByText("Pn3m · 0.91")).toBeInTheDocument();
  });

  it("renders grayed chip for 'dead' state with data-mention-state attribute", () => {
    renderWithProviders(<MentionChip resolved="dead" originalText="Pn3m · 0.91" />);
    const chip = screen.getByText("Pn3m · 0.91");
    expect(chip).toBeInTheDocument();
    expect(chip.closest("[data-mention-state='dead']")).toBeInTheDocument();
  });
});
