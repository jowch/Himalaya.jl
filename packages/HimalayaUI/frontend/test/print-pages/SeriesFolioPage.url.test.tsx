/**
 * FOL P2-3 — the folio controls live in the URL (?q=&filter=&sort=).
 *
 * Separate from SeriesFolioPage.test.tsx because these tests need the REAL
 * react-router (that file mocks useNavigate module-wide): the round-trip pins
 * read the live location, and the replace-semantics pin walks the history
 * stack. Replace is asserted behaviorally: after N keystrokes, ONE back step
 * lands on the entry BELOW /series — keystrokes never became history entries.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation, useNavigate } from "react-router-dom";
import type { SeriesSummary } from "../../src/api";

// ── mock data plane (queries only; the router stays real) ────────────────────
const state = { summaries: [] as SeriesSummary[] };

vi.mock("../../src/queries", () => ({
  useSeriesList: () => ({
    data: state.summaries,
    isLoading: false,
    isError: false,
    isFetching: false,
    refetch: vi.fn(),
  }),
  useSeries: () => ({ data: undefined, isLoading: false }),
  useSeriesTraces: () => ({ data: {}, isLoading: false }),
}));

vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

import { SeriesFolioPage } from "../../src/print/pages/SeriesFolioPage";

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

function LocationSpy(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="loc">{loc.pathname + loc.search}</div>;
}

function BackProbe(): JSX.Element {
  const navigate = useNavigate();
  return (
    <button data-testid="back-probe" onClick={() => navigate(-1)}>
      back
    </button>
  );
}

function renderAt(search: string) {
  return render(
    <MemoryRouter initialEntries={["/elsewhere", `/series${search}`]} initialIndex={1}>
      <Routes>
        <Route
          path="/series"
          element={
            <>
              <SeriesFolioPage />
              <LocationSpy />
              <BackProbe />
            </>
          }
        />
        <Route path="/elsewhere" element={<div data-testid="elsewhere" />} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  state.summaries = [
    summary({ id: 1, title: "LL37 titration lipid 1-2", member_phase_count: 2 }),
    summary({ id: 2, title: "Lipid baselines", member_phase_count: 1 }),
    summary({ id: 3, title: "April vs July cross-exp", member_phase_count: 1, spans_experiments: true }),
  ];
});

describe("SeriesFolioPage URL round-trip (FOL P2-3)", () => {
  it("mounting with ?q=LL37&filter=transition reflects both in the controls", () => {
    renderAt("?q=LL37&filter=transition");
    expect(screen.getByRole("textbox")).toHaveValue("LL37");
    const chips = screen.getAllByTestId("filter-chip");
    const transition = chips.find((c) => c.textContent === "Has transition");
    expect(transition).toHaveAttribute("aria-pressed", "true");
    // Both narrow: only series 1 matches "LL37" AND has a transition.
    expect(screen.getAllByTestId("series-card")).toHaveLength(1);
  });

  it("mounting with ?sort=size selects the Largest segment", () => {
    renderAt("?sort=size");
    expect(screen.getByRole("button", { name: "Largest" })).toHaveAttribute(
      "aria-pressed",
      "true",
    );
  });

  it("typing writes ?q= to the URL", () => {
    renderAt("");
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "Lipid" } });
    expect(screen.getByTestId("loc")).toHaveTextContent("/series?q=Lipid");
  });

  it("keystrokes REPLACE the history entry — one back step leaves /series entirely", () => {
    renderAt("");
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "L" } });
    fireEvent.change(input, { target: { value: "LL" } });
    fireEvent.change(input, { target: { value: "LL37" } });
    expect(screen.getByTestId("loc")).toHaveTextContent("/series?q=LL37");
    // If each keystroke had pushed, back would land on ?q=LL. Replace means
    // the stack is still [/elsewhere, /series] and back exits the folio.
    fireEvent.click(screen.getByTestId("back-probe"));
    expect(screen.getByTestId("elsewhere")).toBeInTheDocument();
  });

  it("clearing the search back to empty removes ?q= (defaults absent from the URL)", () => {
    renderAt("?q=LL37");
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "" } });
    expect(screen.getByTestId("loc")).toHaveTextContent(/^\/series$/);
  });

  it("returning the filter to All removes ?filter=", () => {
    renderAt("?filter=transition");
    const chips = screen.getAllByTestId("filter-chip");
    fireEvent.click(chips.find((c) => c.textContent === "All")!);
    expect(screen.getByTestId("loc")).toHaveTextContent(/^\/series$/);
  });

  it("the SearchInput inline clear also clears ?q=", () => {
    renderAt("?q=LL37");
    fireEvent.click(screen.getByRole("button", { name: "Clear search" }));
    expect(screen.getByTestId("loc")).toHaveTextContent(/^\/series$/);
    expect(screen.getByRole("textbox")).toHaveValue("");
  });

  it("'Show the whole folio' clears q + filter params; the sort param SURVIVES", () => {
    renderAt("?q=zzzznotaseries&filter=transition&sort=size");
    expect(screen.getByText("No series match")).toBeInTheDocument();
    const block = screen.getByTestId("empty-state");
    fireEvent.click(
      within(block).getByRole("button", { name: "Show the whole folio" }),
    );
    // Search + filter reset (params dropped); sort=size keeps its param and
    // its control state — the pinned FOL-EMPTY semantic, now URL-durable.
    expect(screen.getByTestId("loc")).toHaveTextContent(/^\/series\?sort=size$/);
    expect(screen.getAllByTestId("series-card")).toHaveLength(3);
    expect(screen.getByRole("button", { name: "Largest" })).toHaveAttribute(
      "aria-pressed",
      "true",
    );
  });
});
