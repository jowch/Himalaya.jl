import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { SeriesSummary } from "../src/api";
import { SeriesFolioPage } from "../src/pages/SeriesFolioPage";

const h = vi.hoisted(() => ({
  listQ: {} as { data?: SeriesSummary[]; isLoading: boolean; isError: boolean },
}));

vi.mock("../src/queries", () => ({
  useSeriesList: () => h.listQ,
}));

function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 1, title: "Alpha series", description: null, content_hash: "h",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-02 10:00:00",
    author_username: "jc", member_count: 3, member_phases: ["Pn3m"],
    member_phase_count: 1, has_stale_members: false, ...over,
  };
}

function renderAt(path = "/series") {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/series" element={<SeriesFolioPage />} />
          <Route path="/series/:id" element={<div data-testid="builder-marker" />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SeriesFolioPage", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.listQ = { data: undefined, isLoading: false, isError: false };
  });

  it("renders a card per series and a count in the header", () => {
    h.listQ = {
      data: [summary({ id: 1, title: "Alpha series" }), summary({ id: 2, title: "Beta series" })],
      isLoading: false, isError: false,
    };
    renderAt();
    expect(screen.getByTestId("series-folio-page")).toBeInTheDocument();
    expect(screen.getByTestId("series-card-1")).toBeInTheDocument();
    expect(screen.getByTestId("series-card-2")).toBeInTheDocument();
    expect(screen.getByTestId("series-folio-count")).toHaveTextContent("2");
  });

  it("shows the boneyard skeleton while loading", () => {
    h.listQ = { data: undefined, isLoading: true, isError: false };
    const { container } = renderAt();
    // While loading, boneyard mounts its wrapper [data-boneyard=<name>] and
    // hides the real children inside [data-boneyard-content] (it overlays
    // measured bones in a real browser; JSDOM has no layout so it emits the
    // wrapper but zero bone divs). Assert the wrapper is present and that no
    // real card has rendered — the page is in its cold-fetch state.
    expect(container.querySelector('[data-boneyard="series-folio"]')).not.toBeNull();
    expect(screen.queryByTestId("series-card-1")).not.toBeInTheDocument();
  });

  it("shows an empty state when there are zero series", () => {
    h.listQ = { data: [], isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-folio-empty")).toBeInTheDocument();
  });

  it("shows an error state when the query errors", () => {
    h.listQ = { data: undefined, isLoading: false, isError: true };
    renderAt();
    expect(screen.getByTestId("series-folio-error")).toBeInTheDocument();
  });

  it("filters cards by the search box (client-side, title match)", () => {
    h.listQ = {
      data: [summary({ id: 1, title: "Alpha series" }), summary({ id: 2, title: "Beta series" })],
      isLoading: false, isError: false,
    };
    renderAt();
    fireEvent.change(screen.getByTestId("series-folio-search"), { target: { value: "beta" } });
    expect(screen.queryByTestId("series-card-1")).not.toBeInTheDocument();
    expect(screen.getByTestId("series-card-2")).toBeInTheDocument();
  });

  it("navigates to /series/:id when a card is opened", () => {
    h.listQ = { data: [summary({ id: 1 })], isLoading: false, isError: false };
    renderAt();
    fireEvent.click(screen.getByTestId("series-card-1"));
    expect(screen.getByTestId("builder-marker")).toBeInTheDocument();
  });
});
