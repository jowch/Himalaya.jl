import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TracePlot, type TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 1, 1] },
  peaks: [{ id: 1, q: 0.2, source: "auto", intensity: 20 }],
  phase: "Ia3d",
};

describe("TracePlot", () => {
  it("renders axes, a trace line, and peak glyphs for a single trace", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} data-testid="tp" />,
    );
    expect(container.querySelector('svg[data-testid="tp"]')).toBeTruthy();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(1);
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeTruthy();
    expect(container.querySelector('[data-role="axis-left"]')).toBeTruthy();
  });

  it("renders no axes when axes={false} (mini level)", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={120} height={48} axes={false} />,
    );
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeNull();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
  });

  it("calls onAddPeak when clicking empty space", () => {
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlot
        traces={[model]}
        width={500}
        height={300}
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    // Click far from the single peak (q=0.2 ≈ px 274 in plot space). px=148.
    fireEvent.click(container.querySelector("svg")!, {
      clientX: 200,
      clientY: 100,
    });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
  });
});
