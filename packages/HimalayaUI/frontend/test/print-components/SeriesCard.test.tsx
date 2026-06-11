import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SeriesCard } from "../../src/print/components/SeriesCard";
import { TRANSITION, FULL } from "../../src/print/waterfall/waterfall.fixtures";
import type { PhaseSegment } from "../../src/print/ui";

const SEGMENTS: PhaseSegment[] = [
  { phase: "Ia3d" },
  { phase: "Im3m", coexistWith: ["Ia3d"] },
  { phase: "Im3m" },
];

const FULL_SEGMENTS: PhaseSegment[] = [
  { phase: "Ia3d" },
  { phase: "Ia3d" },
  { phase: "Im3m", coexistWith: ["Ia3d"] },
  { phase: "Im3m" },
  { phase: "Lamellar" },
];

const BASE_PROPS = {
  rows: TRANSITION,
  segments: SEGMENTS,
  figLabel: "Fig. 1",
  title: "LL37 titration of lipid 1-2",
  sampleCount: 3,
  variable: "LL37 : lipid ratio",
  provenance: "SSRL · Apr 2026",
  editedLabel: "2 days ago",
  author: "JC",
};

describe("<SeriesCard> renders text content", () => {
  it("renders the title", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText("LL37 titration of lipid 1-2")).toBeInTheDocument();
  });

  it("renders the title as a level-2 heading (folio h1 → card h2, no skip)", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(
      screen.getByRole("heading", { level: 2, name: "LL37 titration of lipid 1-2" }),
    ).toBeInTheDocument();
  });

  it("renders the figLabel", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText("Fig. 1")).toBeInTheDocument();
  });

  it("renders the sample count", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText("3 samples")).toBeInTheDocument();
  });

  it("renders the variable", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText(/LL37 : lipid ratio/)).toBeInTheDocument();
  });

  it("renders '· by {variable}' when a variable is present", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText(/by LL37 : lipid ratio/)).toBeInTheDocument();
  });

  it("DROPS the dangling '· by' clause when variable is empty (item 5)", () => {
    render(<SeriesCard {...BASE_PROPS} variable="" />);
    // No dangling preposition: the meta line is just "{n} samples".
    expect(screen.queryByText(/\bby\b/)).toBeNull();
    expect(screen.getByText("3 samples")).toBeInTheDocument();
  });

  it("DROPS the clause for a whitespace-only variable", () => {
    render(<SeriesCard {...BASE_PROPS} variable="   " />);
    expect(screen.queryByText(/\bby\b/)).toBeNull();
  });
});

describe("<SeriesCard> card-figure", () => {
  it("embeds a card-figure whose data-row-count equals rows.length", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    const fig = screen.getByTestId("card-figure");
    expect(fig).toBeInTheDocument();
    expect(fig).toHaveAttribute("data-row-count", String(TRANSITION.length));
  });

  it("full rows fixture has correct row count", () => {
    render(<SeriesCard {...BASE_PROPS} rows={FULL} segments={FULL_SEGMENTS} />);
    const fig = screen.getByTestId("card-figure");
    expect(fig).toHaveAttribute("data-row-count", String(FULL.length));
  });
});

describe("<SeriesCard> phase strip", () => {
  it("embeds ps-seg cells equal in count to segments.length", () => {
    const { getAllByTestId } = render(<SeriesCard {...BASE_PROPS} />);
    expect(getAllByTestId("ps-seg")).toHaveLength(SEGMENTS.length);
  });

  it("embeds a ps-cap", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByTestId("ps-cap")).toBeInTheDocument();
  });
});

describe("<SeriesCard> notice pill", () => {
  it("notice={tone:'new',count:2} renders a notice-pill with data-tone='new' containing '+2 new match'", () => {
    render(<SeriesCard {...BASE_PROPS} notice={{ tone: "new", count: 2 }} />);
    const pill = screen.getByTestId("notice-pill");
    expect(pill).toHaveAttribute("data-tone", "new");
    expect(pill).toHaveTextContent("+2 new match");
  });

  it("notice={tone:'draft'} renders a notice-pill with data-tone='draft'", () => {
    render(<SeriesCard {...BASE_PROPS} notice={{ tone: "draft" }} />);
    const pill = screen.getByTestId("notice-pill");
    expect(pill).toHaveAttribute("data-tone", "draft");
  });

  it("no notice → no notice-pill rendered", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.queryByTestId("notice-pill")).not.toBeInTheDocument();
  });
});

describe("<SeriesCard> plate lift (FOL P2-1)", () => {
  // DESIGN.md §4: folio cards are Plate Lift surfaces. The lift is the Card
  // primitive's `elevated` variant (the single source of the plate look), so
  // the pin is the variant attribute — semantics, not class strings.
  it("a saved card carries the Card elevated variant (data-elevated='true')", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByTestId("series-card")).toHaveAttribute("data-elevated", "true");
  });

  it("a draft card stays flat — dashed recipe, NO Plate Lift", () => {
    render(<SeriesCard {...BASE_PROPS} draft />);
    expect(screen.getByTestId("series-card")).not.toHaveAttribute("data-elevated");
  });

  it("a clickable card carries the quiet hover door affordance (data-interactive='true')", () => {
    render(<SeriesCard {...BASE_PROPS} onClick={vi.fn()} />);
    expect(screen.getByTestId("series-card")).toHaveAttribute("data-interactive", "true");
  });

  it("a non-clickable card has no door affordance", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByTestId("series-card")).not.toHaveAttribute("data-interactive");
  });
});

describe("<SeriesCard> draft variant", () => {
  it("draft=true → root data-draft='true'", () => {
    render(<SeriesCard {...BASE_PROPS} draft />);
    expect(screen.getByTestId("series-card")).toHaveAttribute("data-draft", "true");
  });

  it("default (no draft prop) → root data-draft='false'", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByTestId("series-card")).toHaveAttribute("data-draft", "false");
  });
});

describe("<SeriesCard> onClick", () => {
  it("fires onClick when the card is clicked", () => {
    const onClick = vi.fn();
    render(<SeriesCard {...BASE_PROPS} onClick={onClick} />);
    fireEvent.click(screen.getByTestId("series-card"));
    expect(onClick).toHaveBeenCalledOnce();
  });

  it("exposes the title as a keyboard-focusable button when onClick is set (FOL-KBD)", () => {
    const onClick = vi.fn();
    render(<SeriesCard {...BASE_PROPS} onClick={onClick} />);
    const button = screen.getByRole("button", {
      name: "LL37 titration of lipid 1-2",
    });
    expect(button).toBeInTheDocument();
  });

  it("clicking the title button fires onClick exactly once (stopPropagation, no double-nav)", () => {
    const onClick = vi.fn();
    render(<SeriesCard {...BASE_PROPS} onClick={onClick} />);
    fireEvent.click(
      screen.getByRole("button", { name: "LL37 titration of lipid 1-2" }),
    );
    expect(onClick).toHaveBeenCalledOnce();
  });

  it("the title button still presents as a level-2 heading (SC-HEAD survives)", () => {
    render(<SeriesCard {...BASE_PROPS} onClick={vi.fn()} />);
    expect(
      screen.getByRole("heading", {
        level: 2,
        name: "LL37 titration of lipid 1-2",
      }),
    ).toBeInTheDocument();
  });

  it("renders the title as plain text (no button) when onClick is absent", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.queryByRole("button")).toBeNull();
    expect(
      screen.getByRole("heading", {
        level: 2,
        name: "LL37 titration of lipid 1-2",
      }),
    ).toBeInTheDocument();
  });
});

describe("<SeriesCard> footer", () => {
  it("renders editedLabel in the footer", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText(/2 days ago/)).toBeInTheDocument();
  });

  it("renders author in the footer", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText(/JC/)).toBeInTheDocument();
  });

  it("renders the provenance string in the footer", () => {
    render(<SeriesCard {...BASE_PROPS} />);
    expect(screen.getByText("SSRL · Apr 2026")).toBeInTheDocument();
  });

  it("renders a ReactNode provenance in the footer", () => {
    render(
      <SeriesCard
        {...BASE_PROPS}
        provenance={<span data-testid="prov-node">April + July</span>}
      />
    );
    expect(screen.getByTestId("prov-node")).toBeInTheDocument();
  });
});
