import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ResidualChart } from "../../src/print/comb/ResidualChart";
import type { CombSeries } from "../../src/print/comb/combModel";

const PN3M: CombSeries = {
  phase: "Pn3m",
  color: "oklch(0.570 0.150 58)",
  rSquared: 0.998,
  teeth: [
    { q: 0.0712, label: "√2", observed: true, residual: 0.010 },  // positive → above baseline
    { q: 0.0872, label: "√3", observed: true, residual: -0.008 }, // negative → below baseline
    { q: 0.1126, label: "√5", observed: false },                  // no residual → no point
    { q: 0.1234, label: "√6", observed: true, residual: 0.041 },  // beyond ±3% → overflow
  ],
};

describe("ResidualChart", () => {
  it("renders one resid-row per assigned/preview series and NO leftover row", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} hovered={PN3M} />);
    expect(container.querySelectorAll('[data-role="resid-row"]').length).toBe(2);
    expect(container.querySelector('[data-row-kind="leftover"]')).toBeNull();
  });

  it("shows the R² in the gutter", () => {
    const { getByText } = render(<ResidualChart assigned={[PN3M]} />);
    expect(getByText(/R².*0\.998/)).toBeTruthy();
  });

  it("draws a tolerance band per row", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    expect(container.querySelectorAll('[data-role="tolerance-band"]').length).toBe(1);
  });

  it("plots a point only for observed reflections with a residual", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    expect(container.querySelectorAll('[data-role="resid-point"]').length).toBe(3); // √2 √3 √6, not √5
  });

  it("places a positive residual above the baseline and a negative one below", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    const baseY = Number(container.querySelector('[data-role="row-baseline"]')!.getAttribute("y1"));
    const pos = container.querySelector('[data-role="resid-point"][data-q="0.0712"] circle')!;
    const neg = container.querySelector('[data-role="resid-point"][data-q="0.0872"] circle')!;
    expect(Number(pos.getAttribute("cy"))).toBeLessThan(baseY);    // above = smaller y
    expect(Number(neg.getAttribute("cy"))).toBeGreaterThan(baseY); // below = larger y
  });

  it("marks an out-of-domain residual as overflow (open dot)", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    const overflow = container.querySelector('[data-role="resid-point"][data-overflow="true"] circle')!;
    expect(overflow.getAttribute("data-q") ?? overflow.parentElement!.getAttribute("data-q")).toBe("0.1234");
    expect(overflow.getAttribute("fill")).toBe("none");
  });

  it("a hovered point goes hot keeping its own colour (no accent)", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} hoveredQ={0.0712} />);
    const hot = container.querySelector('[data-role="resid-point"][data-hot="true"]')!;
    expect(hot.getAttribute("data-q")).toBe("0.0712");
    const circle = hot.querySelector("circle")!;
    expect(circle.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
    expect(circle.getAttribute("stroke")).not.toBe("var(--color-accent)");
  });

  it("fires onHoverQ on point enter/leave", () => {
    const spy = vi.fn();
    const { container } = render(<ResidualChart assigned={[PN3M]} onHoverQ={spy} />);
    const point = container.querySelector('[data-role="resid-point"]')!;
    fireEvent.mouseEnter(point);
    expect(spy).toHaveBeenCalledWith(0.0712);
    fireEvent.mouseLeave(point);
    expect(spy).toHaveBeenCalledWith(undefined);
  });
});
