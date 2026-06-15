import { render, act } from "@testing-library/react";
import { describe, it, expect, vi, afterEach } from "vitest";
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

describe("CardFigure responsive width (FOL-FIGCLIP)", () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  it("fills its container by default (width: 100%), sizing rows from the 340px reference before measurement", () => {
    const { getByTestId, container } = render(<CardFigure rows={TRANSITION} />);
    const root = getByTestId("card-figure") as HTMLElement;
    expect(root.style.width).toBe("100%");
    // Unmeasured (JSDOM's ResizeObserver stub never fires): the 340px mockup
    // reference, minus the 14px+14px pads.
    const row = container.querySelector("[data-row-key]") as HTMLElement;
    expect(row.style.width).toBe("312px");
  });

  it("an explicit width prop pins the root and the rows", () => {
    const { getByTestId, container } = render(<CardFigure rows={TRANSITION} width={300} />);
    expect((getByTestId("card-figure") as HTMLElement).style.width).toBe("300px");
    const row = container.querySelector("[data-row-key]") as HTMLElement;
    expect(row.style.width).toBe("272px"); // 300 − 14 − 14
  });

  it("tracks the measured container width (no clipping at narrow card widths)", () => {
    let cb: ((entries: Array<{ contentRect: { width: number } }>) => void) | undefined;
    class ROStub {
      constructor(c: (entries: Array<{ contentRect: { width: number } }>) => void) {
        cb = c;
      }
      observe(): void {}
      unobserve(): void {}
      disconnect(): void {}
    }
    vi.stubGlobal("ResizeObserver", ROStub);
    const { container } = render(<CardFigure rows={TRANSITION} />);
    expect(cb).toBeDefined();
    // The 1024px-viewport card interior measures ~285px — the figure follows.
    act(() => {
      cb!([{ contentRect: { width: 285 } }]);
    });
    const row = container.querySelector("[data-row-key]") as HTMLElement;
    expect(row.style.width).toBe("257px"); // 285 − 14 − 14: no overflow clip
  });
});
