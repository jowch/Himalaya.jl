import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { PlotFrame } from "../../src/print/plot/PlotFrame";

describe("PlotFrame", () => {
  const margins = { top: 10, right: 10, bottom: 20, left: 30 };

  it("renders an svg and calls render with the computed plot dims", () => {
    const renderSpy = vi.fn(() => null);
    render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        data-testid="frame"
        render={renderSpy}
      />,
    );
    expect(renderSpy).toHaveBeenCalledWith(
      expect.objectContaining({
        width: 400,
        height: 200,
        plotWidth: 360, // 400 - 30 - 10
        plotHeight: 170, // 200 - 10 - 20
        margins,
      }),
    );
  });

  it("renders a testid'd svg element", () => {
    const { container } = render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        data-testid="frame"
        render={() => null}
      />,
    );
    expect(container.querySelector('svg[data-testid="frame"]')).toBeTruthy();
  });

  it("calls onPointerMovePx when pointer moves over the container div", () => {
    const spy = vi.fn();
    const { container } = render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        render={() => null}
        onPointerMovePx={spy}
      />,
    );
    const div = container.querySelector("div")!;
    fireEvent.pointerMove(div, { clientX: 30, clientY: 20 });
    // jsdom getBoundingClientRect returns 0,0; clientX/Y are passed through.
    // We only assert the spy was called (coordinate math is the same as click/wheel).
    expect(spy).toHaveBeenCalled();
  });

  it("names the svg as a graphic when svgRole/svgLabel are given (WCAG 1.1.1)", () => {
    const { container } = render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        data-testid="frame"
        svgRole="img"
        svgLabel="Integration trace: intensity vs q"
        render={() => null}
      />,
    );
    const svg = container.querySelector("svg")!;
    expect(svg.getAttribute("role")).toBe("img");
    expect(svg.getAttribute("aria-label")).toBe("Integration trace: intensity vs q");
  });

  it("leaves the svg nameless when svgRole/svgLabel are omitted (existing consumers byte-identical)", () => {
    const { container } = render(
      <PlotFrame height={200} margins={margins} width={400} render={() => null} />,
    );
    const svg = container.querySelector("svg")!;
    expect(svg.getAttribute("role")).toBeNull();
    expect(svg.getAttribute("aria-label")).toBeNull();
  });

  it("calls onPointerLeave when pointer leaves the container div", () => {
    const leaveSpy = vi.fn();
    const { container } = render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        render={() => null}
        onPointerLeave={leaveSpy}
      />,
    );
    const div = container.querySelector("div")!;
    fireEvent.pointerLeave(div);
    expect(leaveSpy).toHaveBeenCalledTimes(1);
  });
});
