import { render, fireEvent, within } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SeriesPlate } from "../../src/print/components/SeriesPlate";
import { TRANSITION } from "../../src/print/waterfall/waterfall.fixtures";

describe("SeriesPlate", () => {
  it("renders header, waterfall, scale toggle, and the foot", () => {
    const { getByTestId, getByText, getAllByRole } = render(
      <SeriesPlate kicker="Series" title="LL37 titration" subtitle="7 exposures" rows={TRANSITION}
        scale="log" onScaleChange={() => {}} legendPhases={["Pn3m", "Lamellar"]}
        footHint="peaks are light anchors" footNote="offset ×1.0 · log I" />,
    );
    expect(getByTestId("series-plate")).toBeTruthy();
    expect(getByText("LL37 titration")).toBeTruthy();
    expect(getByTestId("waterfall")).toBeTruthy();
    expect(getAllByRole("button").some((b) => b.textContent === "log q")).toBe(true);
    expect(getByText("offset ×1.0 · log I")).toBeTruthy();
  });
  it("does not render a heatmap toggle", () => {
    const { queryByText } = render(
      <SeriesPlate title="x" rows={TRANSITION} scale="log" onScaleChange={() => {}} />,
    );
    expect(queryByText(/heatmap/i)).toBeNull();
  });
  it("renders the actions slot inside the plate (e.g. the figure-export button)", () => {
    const { getByTestId } = render(
      <SeriesPlate title="x" rows={TRANSITION} scale="log" onScaleChange={() => {}}
        actions={<button data-testid="x-action">A</button>} />,
    );
    expect(within(getByTestId("series-plate")).getByTestId("x-action")).toBeInTheDocument();
  });
  it("wires the scale toggle", () => {
    const onScaleChange = vi.fn();
    const { getByText } = render(
      <SeriesPlate title="x" rows={TRANSITION} scale="log" onScaleChange={onScaleChange} />,
    );
    fireEvent.click(getByText("linear q"));
    expect(onScaleChange).toHaveBeenCalledWith("lin");
  });
});
