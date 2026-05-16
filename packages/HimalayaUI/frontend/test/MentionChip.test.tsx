import { describe, it, expect, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionChip } from "../src/components/MentionChip";
import { useAppState } from "../src/state";
import type { ResolvedMention } from "../src/hooks/useMentionResolution";
import type * as api from "../src/api";
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

// ─── Comparison chip (Phase 10) ───────────────────────────────────────────
describe("<MentionChip> — comparison", () => {
  const COMP: api.Comparison = {
    id: 7, title: "DOPE titration",
    description: null,
    // Prefixed `sha256:` form matches the real API payload — this fixture
    // previously stored the bare hex, which masked issue #62 (the drift
    // detector and the picker emit path sliced inconsistently).
    content_hash: "sha256:a1b2c3d4e5f60718",
    created_by: 1,
    created_at: "2026-05-02 10:00:00",
    updated_at: "2026-05-02 10:00:00",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: null,
    members: [
      // Three members — count surfaces in the chip tooltip.
      makeMember(1), makeMember(2), makeMember(3),
    ],
  };

  it("renders title as the chip label", () => {
    const resolved: ResolvedMention = { type: "comparison", data: COMP };
    renderWithProviders(<MentionChip resolved={resolved} />);
    expect(screen.getByText("DOPE titration")).toBeInTheDocument();
  });

  it("exposes data-testid + data-mention-* attrs for E2E selectors", () => {
    const resolved: ResolvedMention = { type: "comparison", data: COMP };
    renderWithProviders(<MentionChip resolved={resolved} tokenHash="a1b2c3d4" />);
    const chip = screen.getByText("DOPE titration").closest("[data-mention-type]")!;
    expect(chip).toHaveAttribute("data-testid", "mention-chip");
    expect(chip).toHaveAttribute("data-mention-type", "comparison");
    expect(chip).toHaveAttribute("data-mention-id", "7");
  });

  it("data-hash-drift='false' when tokenHash matches current content_hash", () => {
    const resolved: ResolvedMention = { type: "comparison", data: COMP };
    renderWithProviders(<MentionChip resolved={resolved} tokenHash="a1b2c3d4" />);
    const chip = screen.getByText("DOPE titration").closest("[data-mention-type]")!;
    expect(chip).toHaveAttribute("data-hash-drift", "false");
  });

  it("data-hash-drift='true' when tokenHash diverges from current content_hash", () => {
    const resolved: ResolvedMention = { type: "comparison", data: COMP };
    renderWithProviders(<MentionChip resolved={resolved} tokenHash="00000000" />);
    const chip = screen.getByText(/DOPE titration/).closest("[data-mention-type]")!;
    expect(chip).toHaveAttribute("data-hash-drift", "true");
    // Visible drift annotation alongside the title.
    expect(screen.getByText(/changed/i)).toBeInTheDocument();
  });

  it("data-hash-drift='false' when no tokenHash is provided", () => {
    // Legacy `[[comparison:7]]` form (no eager hash) — drift indicator off.
    const resolved: ResolvedMention = { type: "comparison", data: COMP };
    renderWithProviders(<MentionChip resolved={resolved} />);
    const chip = screen.getByText("DOPE titration").closest("[data-mention-type]")!;
    expect(chip).toHaveAttribute("data-hash-drift", "false");
  });
});

function makeMember(id: number): api.ComparisonMember {
  return {
    id,
    comparison_id: 7,
    exposure_id: 100 + id,
    display_order: id,
    band_height: 1.0,
    y_offset: 0,
    normalization: "none",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: null,
    is_stale: false,
    created_by: null,
    created_at: null,
  };
}
