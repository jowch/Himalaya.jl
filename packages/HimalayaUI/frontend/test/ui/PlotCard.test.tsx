import { describe, expect, it } from "vitest";
import { render, screen } from "@testing-library/react";
import { PlotCard } from "../../src/ui/PlotCard";

describe("PlotCard", () => {
  it("holds the legend inside the card, above the plot", () => {
    render(
      <PlotCard legend={<span>chips</span>}>
        <svg data-testid="plot" />
      </PlotCard>,
    );
    const legend = screen.getByTestId("plot-card-legend");
    const card = legend.parentElement!;
    expect(legend).toHaveTextContent("chips");
    expect(card).toContainElement(screen.getByTestId("plot"));
    expect(legend.compareDocumentPosition(screen.getByTestId("plot")) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
  });
});
