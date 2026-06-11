import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import type { SeriesSummary, Series, SeriesMember } from "../../src/api";

// ── navigate spy ─────────────────────────────────────────────────────────────
const navigateSpy = vi.fn();
const listRefetch = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return {
    ...actual,
    useNavigate: () => navigateSpy,
  };
});

// ── mock data plane ────────────────────────────────────────────────────────────
const state = {
  summaries: [] as SeriesSummary[],
  seriesById: new Map<number, Series>(),
  loading: false,
  error: false,
  fetching: false,
};

vi.mock("../../src/queries", () => ({
  useSeriesList: () => ({
    data: state.error ? undefined : state.summaries,
    isLoading: state.loading,
    isError: state.error,
    isFetching: state.fetching,
    refetch: listRefetch,
  }),
  useSeries: (id: number | undefined) => ({
    data: id !== undefined ? state.seriesById.get(id) : undefined,
    isLoading: state.loading,
  }),
  useSeriesTraces: (_id: number | undefined) => ({
    data: {},
    isLoading: false,
  }),
}));

// boneyard Skeleton: render children when not loading.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children, loading, fallback }: {
    children: React.ReactNode;
    loading: boolean;
    fallback?: React.ReactNode;
  }) => loading ? <>{fallback}</> : <>{children}</>,
}));

import { SeriesFolioPage } from "../../src/print/pages/SeriesFolioPage";

// ── helpers ────────────────────────────────────────────────────────────────────
function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 1,
    title: "Default series",
    description: null,
    content_hash: "abc123",
    created_by: 1,
    created_at: "2026-01-01T00:00:00Z",
    updated_at: "2026-06-01T00:00:00Z",
    forked_from_id: null,
    forked_at_hash: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: "2026-06-01T00:00:00Z",
    author_username: "JC",
    member_count: 3,
    member_phases: ["Pn3m"],
    member_phase_count: 1,
    has_stale_members: false,
    ordering_variable: "LL37 : lipid ratio",
    spans_experiments: false,
    experiment_name: null,
    ...over,
  };
}

function member(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1,
    series_id: 1,
    exposure_id: 1,
    display_order: 0,
    band_height: 1,
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
    ...over,
  };
}

function seriesDetail(id: number, members: SeriesMember[] = []): Series {
  return {
    id,
    title: "Detail series",
    description: null,
    content_hash: "abc123",
    created_by: 1,
    created_at: "2026-01-01T00:00:00Z",
    updated_at: "2026-06-01T00:00:00Z",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid ratio",
    order_rule: "ascending",
    state: "saved",
    members,
    samples: [],
  };
}

function seed(): void {
  state.summaries = [
    summary({ id: 1, title: "LL37 titration lipid 1-2", member_phase_count: 2, spans_experiments: false, experiment_name: "April 2026 beamtime" }),
    summary({ id: 2, title: "Lipid baselines", member_phase_count: 1, spans_experiments: false }),
    summary({ id: 3, title: "April vs July cross-exp", member_phase_count: 1, spans_experiments: true, content_hash: "" }),
  ];
  state.seriesById = new Map([
    [1, seriesDetail(1, [member({ id: 10, series_id: 1, exposure_id: 1 })])],
    [2, seriesDetail(2, [member({ id: 20, series_id: 2, exposure_id: 2 })])],
    [3, seriesDetail(3, [member({ id: 30, series_id: 3, exposure_id: 3 })])],
  ]);
  state.loading = false;
  state.error = false;
  state.fetching = false;
}

function renderPage() {
  return render(
    <MemoryRouter initialEntries={["/series"]}>
      <SeriesFolioPage />
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  seed();
});

describe("SeriesFolioPage", () => {
  it("renders the FolioHeader with total count and one card per listed series", () => {
    renderPage();
    // Header with testid
    expect(screen.getByTestId("folio-header")).toBeInTheDocument();
    // Shows total count (3)
    expect(screen.getByText("3")).toBeInTheDocument();
    // One card per series
    expect(screen.getAllByTestId("series-card")).toHaveLength(3);
  });

  it("search input narrows the visible cards", () => {
    renderPage();
    // SearchInput's data-testid is on the wrapper div; the actual <input> is inside.
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "LL37" } });
    // Only "LL37 titration lipid 1-2" matches
    expect(screen.getAllByTestId("series-card")).toHaveLength(1);
  });

  it("Has transition filter drops single-phase series", () => {
    renderPage();
    // "Has transition" = member_phase_count > 1; only series 1 qualifies
    const chips = screen.getAllByTestId("filter-chip");
    const transitionChip = chips.find((c) => c.textContent === "Has transition");
    expect(transitionChip).toBeDefined();
    fireEvent.click(transitionChip!);
    // Only series 1 passes (member_phase_count=2); series 2 (1) and 3 (1) drop
    expect(screen.getAllByTestId("series-card")).toHaveLength(1);
  });

  it("clicking a card navigates to /series/:id", () => {
    renderPage();
    const cards = screen.getAllByTestId("series-card");
    fireEvent.click(cards[0]!);
    expect(navigateSpy).toHaveBeenCalledWith("/series/1");
  });

  it("a single-experiment card names its beamtime; a null-provenance card stays silent (FOL P2-2)", () => {
    renderPage();
    const cards = screen.getAllByTestId("series-card");
    // Series 1 carries experiment_name — the footer names the beamtime.
    expect(screen.getByText("April 2026 beamtime")).toBeInTheDocument();
    // Series 2 has experiment_name: null and does not span — no provenance
    // text appears on its card (the footer-left slot is empty).
    const card2 = cards.find((c) => within(c).queryByText("Lipid baselines") !== null);
    expect(card2).toBeDefined();
    expect(within(card2!).queryByText("April 2026 beamtime")).toBeNull();
    expect(within(card2!).queryByText(/cross-experiment/)).toBeNull();
  });

  it("a draft series (content_hash === '') shows the Draft pill", () => {
    renderPage();
    // Series 3 is draft (content_hash="")
    // The NoticePill renders "Draft" text for tone="draft"
    expect(screen.getByText("Draft")).toBeInTheDocument();
  });

  it("shows the empty state when the filter matches nothing", () => {
    renderPage();
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "xxxxxxxxxnotaseriesname" } });
    expect(screen.getByTestId("gallery-empty")).toBeInTheDocument();
    expect(screen.getByTestId("empty-state")).toBeInTheDocument();
  });

  it("a genuinely empty folio shows the first-run state with a door to the contact sheet (FOL-EMPTY-FIRSTRUN)", () => {
    state.summaries = [];
    state.seriesById = new Map();
    renderPage();
    // Honest first-run copy, not the filtered no-match masquerade.
    expect(screen.getByText("No series yet")).toBeInTheDocument();
    expect(screen.queryByText("No series match")).toBeNull();
    expect(
      screen.queryByText("Clear the search or filter to see the whole folio."),
    ).toBeNull();
    // The action is a door to the creation path.
    const block = screen.getByTestId("empty-state");
    fireEvent.click(
      within(block).getByRole("button", { name: "Open the contact sheet" }),
    );
    expect(navigateSpy).toHaveBeenCalledWith("/samples");
  });

  it("zero summaries plus a typed search still shows the first-run state (a search over nothing is still nothing saved)", () => {
    state.summaries = [];
    state.seriesById = new Map();
    renderPage();
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "LL37" } });
    expect(screen.getByText("No series yet")).toBeInTheDocument();
    expect(screen.queryByText("No series match")).toBeNull();
  });

  it("the no-match state offers a one-click clear that shows the whole folio again", () => {
    renderPage();
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "xxxxxxxxxnotaseriesname" } });
    const chips = screen.getAllByTestId("filter-chip");
    const transitionChip = chips.find((c) => c.textContent === "Has transition");
    fireEvent.click(transitionChip!);
    // Dirty the sort too: clearing must NOT touch it.
    fireEvent.click(screen.getByRole("button", { name: "Largest" }));
    expect(screen.getByText("No series match")).toBeInTheDocument();
    const block = screen.getByTestId("empty-state");
    fireEvent.click(
      within(block).getByRole("button", { name: "Show the whole folio" }),
    );
    // Search and filter both reset; all three cards return; sort survives.
    expect(screen.getAllByTestId("series-card")).toHaveLength(3);
    expect(screen.getByRole("textbox")).toHaveValue("");
    expect(
      screen.getByRole("button", { name: "Largest" }),
    ).toHaveAttribute("aria-pressed", "true");
  });

  it("shows an honest error surface when the list fetch fails, not 'No series match'", () => {
    state.error = true;
    renderPage();
    // Distinct error copy — not the zero-results "No series match" line.
    expect(screen.getByText("Couldn't load the folio")).toBeInTheDocument();
    expect(screen.queryByText("No series match")).toBeNull();
    // No cards / no controls Gallery in the error branch.
    expect(screen.queryAllByTestId("series-card")).toHaveLength(0);
  });

  it("the error surface offers a retry control wired to refetch (FOL-ERR)", () => {
    state.error = true;
    renderPage();
    const block = screen.getByTestId("empty-state");
    // The body states what happened; the control embodies the way forward.
    expect(screen.getByText("The series list failed to load.")).toBeInTheDocument();
    expect(screen.queryByText(/try reloading/i)).toBeNull();
    fireEvent.click(within(block).getByRole("button", { name: "Try again" }));
    expect(listRefetch).toHaveBeenCalled();
  });

  it("the retry control is disabled while the refetch is in flight", () => {
    state.error = true;
    state.fetching = true;
    renderPage();
    const block = screen.getByTestId("empty-state");
    expect(within(block).getByRole("button", { name: "Try again" })).toBeDisabled();
  });
});
