import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, within, act } from "@testing-library/react";
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
  detailError: false,
  tracesError: false,
};

// Per-card hook spies (FOL-N+1): record the id each mount passes so tests can
// pin that off-viewport cards pass `undefined` (the enabled:false gate — no
// detail/trace fetch) until they near the viewport.
const { seriesSpy, tracesSpy, detailRefetch, tracesRefetch } = vi.hoisted(() => ({
  seriesSpy: vi.fn<(id: number | undefined) => void>(),
  tracesSpy: vi.fn<(id: number | undefined) => void>(),
  detailRefetch: vi.fn(),
  tracesRefetch: vi.fn(),
}));

vi.mock("../../src/queries", () => ({
  useSeriesList: () => ({
    data: state.error ? undefined : state.summaries,
    isLoading: state.loading,
    isError: state.error,
    isFetching: state.fetching,
    refetch: listRefetch,
  }),
  useSeries: (id: number | undefined) => {
    seriesSpy(id);
    return {
      data: id !== undefined && !state.detailError ? state.seriesById.get(id) : undefined,
      isLoading: state.loading,
      isError: state.detailError,
      refetch: detailRefetch,
    };
  },
  useSeriesTraces: (id: number | undefined) => {
    tracesSpy(id);
    return {
      data: state.tracesError ? undefined : {},
      isLoading: false,
      isError: state.tracesError,
      refetch: tracesRefetch,
    };
  },
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
  state.detailError = false;
  state.tracesError = false;
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

afterEach(() => {
  vi.unstubAllGlobals();
});

/** Find the rendered card containing `title`. */
function cardWithTitle(title: string): HTMLElement {
  const card = screen
    .getAllByTestId("series-card")
    .find((c) => within(c).queryByText(title) !== null);
  expect(card).toBeDefined();
  return card!;
}

/** Controllable IntersectionObserver stub: nothing intersects until a test
 *  fires an instance's callback. Instance order == card render order. */
class IOStub {
  static instances: IOStub[] = [];
  callback: (entries: Array<{ isIntersecting: boolean }>) => void;
  constructor(cb: (entries: Array<{ isIntersecting: boolean }>) => void) {
    this.callback = cb;
    IOStub.instances.push(this);
  }
  observe(): void {}
  unobserve(): void {}
  disconnect(): void {}
}

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

  it("shows a visible 'Sort' label left of the sort control, aria-hidden so SR isn't double-announced (FOL-SORT)", () => {
    renderPage();
    const label = screen.getByText("Sort");
    expect(label).toBeInTheDocument();
    expect(label).toHaveAttribute("aria-hidden", "true");
    // The control still names itself for screen readers.
    expect(screen.getByLabelText("Sort series")).toBeInTheDocument();
  });

  it("search input narrows the visible cards", () => {
    renderPage();
    // SearchInput's data-testid is on the wrapper div; the actual <input> is inside.
    const input = screen.getByRole("textbox");
    // "titration" appears only in series 1's title (the seed's default
    // ordering_variable "LL37 : lipid ratio" carries "LL37" across all three,
    // so a title-unique token is the clean single-match probe under the
    // broadened search haystack — F1).
    fireEvent.change(input, { target: { value: "titration" } });
    expect(screen.getAllByTestId("series-card")).toHaveLength(1);
  });

  it("the search field carries a real accessible name, not the disappearing placeholder (FOL-MISC / WCAG 3.3.2)", () => {
    renderPage();
    // The accname must come from an explicit label, NOT the placeholder
    // "Search series…" (which vanishes once the user types). Assert the exact
    // aria-label string with no ellipsis — passes only when a real label wins
    // over the placeholder fallback.
    expect(
      screen.getByRole("textbox", { name: "Search series" }),
    ).toBeInTheDocument();
  });

  it("broadened search matches a member-phase token, not just the title (F1)", () => {
    state.summaries = [
      summary({ id: 1, title: "Cubic run", member_phases: ["Pn3m", "Im3m"] }),
      summary({ id: 2, title: "Flat run", member_phases: ["Lamellar"] }),
    ];
    state.seriesById = new Map([
      [1, seriesDetail(1, [member({ id: 10, series_id: 1, exposure_id: 1 })])],
      [2, seriesDetail(2, [member({ id: 20, series_id: 2, exposure_id: 2 })])],
    ]);
    renderPage();
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "Im3m" } });
    expect(screen.getAllByTestId("series-card")).toHaveLength(1);
    expect(screen.getByText("Cubic run")).toBeInTheDocument();
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

describe("SeriesFolioPage flexibility pass (F2/F3/F4)", () => {
  it("renders a live result count: 'N series' at rest", () => {
    renderPage();
    expect(screen.getByTestId("folio-result-count")).toHaveTextContent("3 series");
  });

  it("the live count switches to 'Showing N of M' once a search narrows the wall", () => {
    renderPage();
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "titration" } });
    expect(screen.getByTestId("folio-result-count")).toHaveTextContent("Showing 1 of 3");
  });

  it("the sort-direction toggle flips the active sort (size: largest ⇄ smallest first)", () => {
    state.summaries = [
      summary({ id: 1, title: "Small", member_count: 2 }),
      summary({ id: 2, title: "Big", member_count: 9 }),
    ];
    state.seriesById = new Map([
      [1, seriesDetail(1, [member({ id: 10, series_id: 1, exposure_id: 1 })])],
      [2, seriesDetail(2, [member({ id: 20, series_id: 2, exposure_id: 2 })])],
    ]);
    renderPage();
    fireEvent.click(screen.getByRole("button", { name: "Largest" }));
    // Default (desc): Big first.
    let cards = screen.getAllByTestId("series-card");
    expect(within(cards[0]!).getByText("Big")).toBeInTheDocument();
    // Flip to ascending → Small first.
    fireEvent.click(screen.getByTestId("folio-sort-dir"));
    cards = screen.getAllByTestId("series-card");
    expect(within(cards[0]!).getByText("Small")).toBeInTheDocument();
  });

  it("the direction toggle reflects state via aria-pressed (pressed == ascending)", () => {
    renderPage();
    // Default sort 'recent' → desc; toggle not pressed.
    const toggle = screen.getByTestId("folio-sort-dir");
    expect(toggle).toHaveAttribute("aria-pressed", "false");
    fireEvent.click(toggle);
    expect(screen.getByTestId("folio-sort-dir")).toHaveAttribute("aria-pressed", "true");
  });

  it("renders the new-series ghost tile at the end of a non-empty wall, linking to /series/new", () => {
    renderPage();
    const tile = screen.getByTestId("new-series-tile");
    expect(tile).toBeInTheDocument();
    fireEvent.click(tile);
    expect(navigateSpy).toHaveBeenCalledWith("/series/new");
  });

  it("the ghost tile is NOT shown when no cards match (the empty state wins)", () => {
    renderPage();
    fireEvent.change(screen.getByRole("textbox"), {
      target: { value: "xxxxxxxxxnotaseriesname" },
    });
    expect(screen.queryByTestId("new-series-tile")).toBeNull();
    expect(screen.getByTestId("gallery-empty")).toBeInTheDocument();
  });

  it("a figure-error card offers a per-card retry that re-runs the figure fetch (F4)", () => {
    state.detailError = true;
    renderPage();
    // JSDOM has no IntersectionObserver → cards are near → error tiles render.
    const retry = screen.getAllByRole("button", { name: "Try again" })[0]!;
    fireEvent.click(retry);
    expect(detailRefetch).toHaveBeenCalled();
  });

  it("the filter chips carry an explanatory tooltip (F5)", () => {
    renderPage();
    const chips = screen.getAllByTestId("filter-chip");
    const all = chips.find((c) => c.textContent === "All")!;
    const transition = chips.find((c) => c.textContent === "Has transition")!;
    const cross = chips.find((c) => c.textContent === "Cross-experiment")!;
    // FOL-ALLCHIP-TOOLTIP: the "All" chip is no longer the one chip without an
    // explanation — every filter chip carries a tooltip.
    expect(all).toHaveAttribute("title", "Every saved series, no filter applied");
    expect(transition).toHaveAttribute("title", "Series whose members span more than one phase");
    expect(cross).toHaveAttribute("title", "Series whose members span more than one experiment");
  });

  it("the filter-chip explanation is also reachable by keyboard/SR via aria-describedby (FIX 2)", () => {
    renderPage();
    // Accessible NAME is the label only (description must not bleed in).
    const all = screen.getByRole("button", { name: "All" });
    const transition = screen.getByRole("button", { name: "Has transition" });
    const cross = screen.getByRole("button", { name: "Cross-experiment" });
    const aDesc = document.getElementById(all.getAttribute("aria-describedby")!);
    const tDesc = document.getElementById(transition.getAttribute("aria-describedby")!);
    const cDesc = document.getElementById(cross.getAttribute("aria-describedby")!);
    expect(aDesc).toHaveTextContent("Every saved series, no filter applied");
    expect(tDesc).toHaveTextContent("Series whose members span more than one phase");
    expect(cDesc).toHaveTextContent("Series whose members span more than one experiment");
    expect(aDesc).toHaveClass("sr-only");
    expect(tDesc).toHaveClass("sr-only");
    expect(cDesc).toHaveClass("sr-only");
  });
});

describe("SeriesFolioPage stable fig numbers (FOL-FIGNUM)", () => {
  it("filtering to a single card keeps its corpus-stable number", () => {
    renderPage();
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "Lipid baselines" } });
    const cards = screen.getAllByTestId("series-card");
    expect(cards).toHaveLength(1);
    // Series 2 is the second committed series by id — it stays "Fig. 2" even
    // when it is the only card on screen (never renumbered to "Fig. 1").
    expect(within(cards[0]!).getByText("Fig. 2")).toBeInTheDocument();
    expect(screen.queryByText("Fig. 1")).toBeNull();
  });

  it("changing the sort does not renumber the cards", () => {
    // Make series 2 the largest so "Largest" reorders it to the front.
    state.summaries[1] = summary({ id: 2, title: "Lipid baselines", member_count: 9 });
    renderPage();
    fireEvent.click(screen.getByRole("button", { name: "Largest" }));
    const first = screen.getAllByTestId("series-card")[0]!;
    expect(within(first).getByText("Lipid baselines")).toBeInTheDocument();
    expect(within(first).getByText("Fig. 2")).toBeInTheDocument();
    expect(within(cardWithTitle("LL37 titration lipid 1-2")).getByText("Fig. 1")).toBeInTheDocument();
  });

  it("a draft consumes no fig number — committed series after it stay densely numbered", () => {
    // Seed: ids 1, 2 committed; 3 draft. Add committed id 4 — it must be
    // "Fig. 3" (the draft never held a number), not "Fig. 4".
    state.summaries.push(summary({ id: 4, title: "Fourth series" }));
    state.seriesById.set(4, seriesDetail(4, [member({ id: 40, series_id: 4, exposure_id: 4 })]));
    renderPage();
    expect(within(cardWithTitle("Fourth series")).getByText("Fig. 3")).toBeInTheDocument();
  });

  it("the draft card itself shows no fig number", () => {
    renderPage();
    expect(within(cardWithTitle("April vs July cross-exp")).queryByText(/^Fig\./)).toBeNull();
  });
});

describe("SeriesFolioPage viewport-lazy card data (FOL-N+1)", () => {
  it("off-viewport cards do not fetch detail/traces; nearing the viewport starts the fetch", () => {
    vi.stubGlobal("IntersectionObserver", IOStub);
    IOStub.instances = [];
    renderPage();
    // All card chrome renders from the LIST summary alone…
    expect(screen.getAllByTestId("series-card")).toHaveLength(3);
    // …and no card has mounted an enabled detail/trace query yet.
    expect(seriesSpy.mock.calls.every(([id]) => id === undefined)).toBe(true);
    expect(tracesSpy.mock.calls.every(([id]) => id === undefined)).toBe(true);
    expect(IOStub.instances.length).toBe(3);

    // First card nears the viewport → its queries enable with its id.
    act(() => {
      IOStub.instances[0]!.callback([{ isIntersecting: true }]);
    });
    expect(seriesSpy.mock.calls.some(([id]) => id === 1)).toBe(true);
    expect(tracesSpy.mock.calls.some(([id]) => id === 1)).toBe(true);
    // The other two cards still have not fetched.
    expect(seriesSpy.mock.calls.some(([id]) => id === 2 || id === 3)).toBe(false);
    expect(tracesSpy.mock.calls.some(([id]) => id === 2 || id === 3)).toBe(false);
  });

  it("without IntersectionObserver support, every card fetches immediately (fail-open)", () => {
    // JSDOM has no IntersectionObserver — this is the default environment.
    renderPage();
    for (const id of [1, 2, 3]) {
      expect(tracesSpy.mock.calls.some(([got]) => got === id)).toBe(true);
      expect(seriesSpy.mock.calls.some(([got]) => got === id)).toBe(true);
    }
  });
});

describe("SeriesFolioPage loading skeleton (FOL-BONES)", () => {
  it("the loading state renders the skeleton card grid, not bare text", () => {
    state.loading = true;
    renderPage();
    // The mocked Skeleton renders `fallback` while loading: the house standard
    // is a card-shaped placeholder grid, never a bare "Loading series…" line.
    expect(screen.getByTestId("folio-bones-fallback")).toBeInTheDocument();
    expect(screen.queryByText("Loading series…")).toBeNull();
    expect(screen.queryAllByTestId("series-card")).toHaveLength(0);
  });
});

describe("SeriesFolioPage FOL-HONEST-DERIVED: figureState wiring", () => {
  it("detail-error on a near card → card shows card-figure-error / honest copy, NOT 'No clear phase'", () => {
    // JSDOM has no IntersectionObserver → all cards are near immediately (fail-open).
    state.detailError = true;
    renderPage();
    // All cards should now show the error figure state
    const cards = screen.getAllByTestId("series-card");
    expect(cards.length).toBeGreaterThan(0);
    // card-figure-error appears on every card since detail is errored for all
    expect(screen.getAllByTestId("card-figure-error").length).toBeGreaterThan(0);
    // Each card shows honest text — use getAllByText since 3 cards all show it
    expect(screen.getAllByText("Couldn't load this figure").length).toBeGreaterThan(0);
    // The false "No clear phase" caption must NOT appear anywhere
    expect(screen.queryByText("No clear phase")).not.toBeInTheDocument();
  });

  it("detail-error card still renders chrome (title, fig label, meta) from list summary", () => {
    state.detailError = true;
    renderPage();
    // Chrome is derived from the healthy list summary and must persist
    expect(screen.getByText("LL37 titration lipid 1-2")).toBeInTheDocument();
    expect(screen.getByText("Fig. 1")).toBeInTheDocument();
  });

  it("traces-error on a near card → card shows card-figure-error / honest copy, NOT 'No clear phase'", () => {
    state.tracesError = true;
    renderPage();
    expect(screen.getAllByTestId("card-figure-error").length).toBeGreaterThan(0);
    // Each card shows honest text — multiple cards, use getAllByText
    expect(screen.getAllByText("Couldn't load this figure").length).toBeGreaterThan(0);
    expect(screen.queryByText("No clear phase")).not.toBeInTheDocument();
  });

  it("header shows '—' (not '0') while the list is loading", () => {
    state.loading = true;
    renderPage();
    // The Skeleton renders fallback when loading, but the header is OUTSIDE the
    // Skeleton component — it always renders. Count must show '—', not '0'.
    const countEl = screen.getByTestId("folio-count");
    expect(countEl).toHaveTextContent("—");
    expect(countEl).not.toHaveTextContent("0");
  });

  it("header shows '—' (not '0') in the list-error branch", () => {
    state.error = true;
    renderPage();
    // The error branch renders FolioHeader — it must use null, not 0.
    const countEl = screen.getByTestId("folio-count");
    expect(countEl).toHaveTextContent("—");
    expect(countEl).not.toHaveTextContent("0");
  });

  it("normal load with 3 series: count still shows '3' (regression guard)", () => {
    renderPage();
    expect(screen.getByTestId("folio-count")).toHaveTextContent("3");
  });
});
