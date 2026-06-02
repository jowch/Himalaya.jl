import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CardFigure } from "../../src/print/waterfall/CardFigure";
import { TRANSITION } from "../../src/print/waterfall/waterfall.fixtures";
import { phaseColor } from "../../src/phases";

/** Read a color back through a DOM element so the assertion compares the
 *  JSDOM-canonicalized form (oklch(0.570 …) → oklch(0.57 …)), never a literal. */
function readBackColor(color: string): string {
  const span = document.createElement("span");
  span.style.color = color;
  return span.style.color;
}

describe("CardFigure", () => {
  it("exposes the row count on the root", () => {
    const { getByTestId } = render(<CardFigure rows={TRANSITION} />);
    expect(getByTestId("card-figure").getAttribute("data-row-count")).toBe("3");
  });

  it("renders exactly one element per row, with the row's phase", () => {
    const { container } = render(<CardFigure rows={TRANSITION} />);
    const rowEls = Array.from(container.querySelectorAll("[data-row-key]"));
    expect(rowEls.length).toBe(3);
    const phases = rowEls.map((el) => el.getAttribute("data-phase"));
    expect(phases).toEqual(TRANSITION.map((r) => r.phase ?? "none"));
  });

  it("carries the phase color on the row element (round-tripped, not literal)", () => {
    const { container } = render(<CardFigure rows={TRANSITION} />);
    const firstPhase = TRANSITION[0]!.phase!;
    const rowEls = Array.from(container.querySelectorAll("[data-row-key]")) as HTMLElement[];
    const firstRow = rowEls[0]!;
    expect(firstRow.style.color).toBe(readBackColor(phaseColor(firstPhase)));
  });
});
