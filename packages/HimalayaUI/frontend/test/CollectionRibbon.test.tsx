import { describe, expect, it } from "vitest";
import { render, screen } from "@testing-library/react";
import { CollectionRibbon } from "../src/ui/CollectionRibbon";
import type { CollectionTimeline } from "../src/api";

const TL: CollectionTimeline = {
  start: "2026-04-25T10:00:00",
  samples: [[0, 300, 1], [420, 900, 0], [4500, 5400, 0]],
  beam_off: [5400],
};

describe("CollectionRibbon", () => {
  it("draws one stretch per sample run, reviewed ones marked", () => {
    const { container } = render(<CollectionRibbon timeline={TL} review={{ reviewed: 1, total: 4 }} />);
    const spans = Array.from(container.querySelectorAll("rect[data-reviewed]"));
    expect(spans.map((r) => r.getAttribute("data-reviewed"))).toEqual(["true", "false", "false"]);
    expect(container.querySelectorAll("rect[data-beam-off]")).toHaveLength(1);
  });

  it("places a stretch in proportion to when it was measured", () => {
    const { container } = render(<CollectionRibbon timeline={TL} review={{ reviewed: 1, total: 4 }} />);
    const last = container.querySelectorAll("rect[data-reviewed]")[2]!;
    expect(Number(last.getAttribute("x"))).toBeCloseTo(1000 * 4500 / 5400, 3);
  });

  it("captions the start, the length and the review count of every sample, timed or not", () => {
    render(<CollectionRibbon timeline={TL} review={{ reviewed: 1, total: 4 }} />);
    expect(screen.getByText("Sat 25 Apr 10:00 · 1 h 30 min")).toBeInTheDocument();
    expect(screen.getByText("3 to review")).toBeInTheDocument();
    expect(screen.getByRole("img")).toHaveAccessibleName(/1 of 4 samples reviewed/);
  });

  it("keeps its row and the review count with an empty track when no image has a time", () => {
    render(<CollectionRibbon timeline={null} review={{ reviewed: 0, total: 2 }} />);
    expect(screen.getByTestId("collection-ribbon")).toHaveAttribute("data-empty", "true");
    expect(screen.getByText("No image times yet")).toBeInTheDocument();
    expect(screen.getByText("2 to review")).toBeInTheDocument();
  });

  it("says nothing about review when no sample can be reviewed", () => {
    render(<CollectionRibbon timeline={TL} review={{ reviewed: 0, total: 0 }} />);
    expect(screen.queryByText("All reviewed")).toBeNull();
    expect(screen.queryByText(/to review/)).toBeNull();
    expect(screen.getByRole("img")).not.toHaveAccessibleName(/reviewed/);
  });

  it("draws a span with no review state apart from reviewed and unreviewed ones", () => {
    const tl: CollectionTimeline = { start: "2026-04-25T10:00:00", samples: [[0, 300, 2], [420, 900, 1], [4500, 5400, 0]], beam_off: [] };
    const { container } = render(<CollectionRibbon timeline={tl} review={{ reviewed: 1, total: 2 }} />);
    const spans = Array.from(container.querySelectorAll("rect[data-reviewed]"));
    expect(spans.map((r) => r.getAttribute("data-reviewable"))).toEqual(["false", "true", "true"]);
    expect(spans.map((r) => r.getAttribute("data-reviewed"))).toEqual(["false", "true", "false"]);
  });
});
