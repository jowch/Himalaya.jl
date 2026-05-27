import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import type { SeriesSummary } from "../src/api";
import { SeriesFolioCard } from "../src/components/SeriesFolioCard";

function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 4, title: "Lipid dose response", description: null,
    content_hash: "h", created_by: 1, created_at: null,
    updated_at: null, forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-02 10:00:00",
    author_username: "jc", member_count: 3,
    member_phases: ["Pn3m", "Lamellar"], member_phase_count: 2,
    has_stale_members: false, ...over,
  };
}

function renderCard(s: SeriesSummary, onOpen = vi.fn()) {
  render(
    <MemoryRouter>
      <SeriesFolioCard series={s} onOpen={onOpen} />
    </MemoryRouter>,
  );
  return { onOpen };
}

describe("SeriesFolioCard", () => {
  it("renders the title, author and member count", () => {
    renderCard(summary());
    expect(screen.getByText("Lipid dose response")).toBeInTheDocument();
    expect(screen.getByTestId("series-card-4")).toHaveAttribute(
      "data-member-count", "3",
    );
    expect(screen.getByText(/jc/)).toBeInTheDocument();
  });

  it("marks a draft (content_hash === '') via data-draft", () => {
    renderCard(summary({ content_hash: "" }));
    expect(screen.getByTestId("series-card-4")).toHaveAttribute("data-draft", "true");
  });

  it("marks stale members via data-stale", () => {
    renderCard(summary({ has_stale_members: true }));
    expect(screen.getByTestId("series-card-4")).toHaveAttribute("data-stale", "true");
  });

  it("renders one swatch per top phase with a +N more when overflowing", () => {
    renderCard(summary({ member_phases: ["Pn3m", "Lamellar", "Im3m"], member_phase_count: 5 }));
    expect(screen.getAllByTestId(/series-card-4-swatch-/)).toHaveLength(3);
    expect(screen.getByText("+2 more")).toBeInTheDocument();
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
