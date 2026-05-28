import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import type { SeriesSummary, Series, SeriesMember } from "../src/api";
import { SeriesFolioCard } from "../src/components/SeriesFolioCard";

// The card pulls full per-sample detail (ordered members + snapshots) to draw
// the live miniature + phase strip — mock that hook.
const h = vi.hoisted(() => ({
  detailQ: {} as { data?: Series; isLoading: boolean },
}));
vi.mock("../src/queries", () => ({
  useSeries: () => h.detailQ,
}));

function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 4, title: "Lipid dose response", description: null,
    content_hash: "h", created_by: 1, created_at: null,
    updated_at: "2026-05-25 10:00:00", forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-25 10:00:00",
    author_username: "jc", member_count: 3,
    member_phases: ["Pn3m", "Lamellar"], member_phase_count: 2,
    has_stale_members: false, ...over,
  };
}

function member(phase: string | null, order: number): SeriesMember {
  return {
    id: order + 1, series_id: 4, exposure_id: order + 100, display_order: order,
    band_height: 1, y_offset: 0, normalization: "max", color_override: null,
    label_override: null, q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase === null ? null : {
      effective_peaks: [{ id: 1, q: 0.04, intensity: 1, sharpness: 1, source: "auto" as const }],
      confirmed_index: { id: order + 1, phase, lattice_d: 100, r_squared: 0.99, ngc: 3, peak_ids: [1] },
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: 1, created_at: null,
  };
}

function detail(over: Partial<Series> = {}): Series {
  return {
    id: 4, title: "Lipid dose response", description: null, content_hash: "h",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid ratio", order_rule: "ascending", state: "ready",
    members: [member("Pn3m", 0), member("Lamellar", 1)],
    samples: [], ...over,
  };
}

const NOW = new Date("2026-05-27T10:00:00Z");

function renderCard(s: SeriesSummary, props: Partial<{ onOpen: () => void; figNumber: number; newMatches: number }> = {}) {
  const onOpen = props.onOpen ?? vi.fn();
  render(
    <MemoryRouter>
      <SeriesFolioCard
        series={s}
        onOpen={onOpen}
        now={NOW}
        figNumber={props.figNumber}
        newMatches={props.newMatches}
      />
    </MemoryRouter>,
  );
  return { onOpen };
}

describe("SeriesFolioCard", () => {
  beforeEach(() => {
    h.detailQ = { data: undefined, isLoading: false };
  });

  it("renders the title, author and member count", () => {
    renderCard(summary());
    expect(screen.getByText("Lipid dose response")).toBeInTheDocument();
    expect(screen.getByTestId("series-card-4")).toHaveAttribute("data-member-count", "3");
    expect(screen.getByText(/jc/)).toBeInTheDocument();
  });

  it("falls back to the phase-swatch strip while detail is loading", () => {
    h.detailQ = { data: undefined, isLoading: true };
    renderCard(summary());
    expect(screen.getByTestId("series-card-4-swatches")).toBeInTheDocument();
    expect(screen.queryByTestId("series-mini-waterfall")).not.toBeInTheDocument();
  });

  it("renders the live miniature + phase strip once detail loads", () => {
    h.detailQ = { data: detail(), isLoading: false };
    renderCard(summary());
    expect(screen.getByTestId("series-mini-waterfall")).toBeInTheDocument();
    expect(screen.getByTestId("series-phase-strip")).toBeInTheDocument();
    expect(screen.queryByTestId("series-card-4-swatches")).not.toBeInTheDocument();
  });

  it("shows a Draft pill for an uncommitted series", () => {
    renderCard(summary({ content_hash: "" }));
    const card = screen.getByTestId("series-card-4");
    expect(card).toHaveAttribute("data-draft", "true");
    expect(screen.getByTestId("series-card-4-pill")).toHaveTextContent(/draft/i);
  });

  it("shows a '+N new match' pill when newMatches > 0 on a committed series", () => {
    renderCard(summary(), { newMatches: 2 });
    expect(screen.getByTestId("series-card-4-pill")).toHaveTextContent("+2 new match");
  });

  it("shows the figure number kicker on a committed series", () => {
    renderCard(summary(), { figNumber: 3 });
    expect(screen.getByTestId("series-card-4-kicker")).toHaveTextContent("Fig. 3");
  });

  it("shows a Recipe kicker on a draft", () => {
    renderCard(summary({ content_hash: "" }), { figNumber: 3 });
    expect(screen.getByTestId("series-card-4-kicker")).toHaveTextContent(/recipe/i);
  });

  it("renders a footer with the author and an edited timestamp", () => {
    renderCard(summary({ updated_at: "2026-05-25 10:00:00" }));
    const foot = screen.getByTestId("series-card-4-foot");
    expect(foot).toHaveTextContent("jc");
    expect(foot).toHaveTextContent(/2 days ago/i);
  });

  it("surfaces the ordering variable in the meta line when detail provides one", () => {
    h.detailQ = { data: detail({ ordering_variable: "LL37 : lipid ratio" }), isLoading: false };
    renderCard(summary());
    expect(screen.getByTestId("series-card-4-meta")).toHaveTextContent("LL37 : lipid ratio");
  });

  it("marks stale members via data-stale", () => {
    renderCard(summary({ has_stale_members: true }));
    expect(screen.getByTestId("series-card-4")).toHaveAttribute("data-stale", "true");
  });

  it("calls onOpen with the id when clicked", () => {
    const { onOpen } = renderCard(summary());
    fireEvent.click(screen.getByTestId("series-card-4"));
    expect(onOpen).toHaveBeenCalledWith(4);
  });

  it("shows an untitled fallback when title is empty", () => {
    renderCard(summary({ title: "" }));
    expect(screen.getByText(/untitled series/i)).toBeInTheDocument();
  });
});
