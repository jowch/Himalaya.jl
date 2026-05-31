import { render } from "@testing-library/react";
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
});
