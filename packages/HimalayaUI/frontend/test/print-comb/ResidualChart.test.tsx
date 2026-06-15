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

  it("shows the R² in the gutter, prefixed so it reads as fit quality", () => {
    const { getByText } = render(<ResidualChart assigned={[PN3M]} />);
    expect(getByText(/fit R² 0\.998/)).toBeTruthy();
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

  it("draws a HOLLOW circle at its actual value for an out-of-tolerance but on-scale residual (band < |res| ≤ domain)", () => {
    const s: CombSeries = {
      phase: "Pn3m", color: "oklch(0.570 0.150 58)",
      teeth: [{ q: 0.09, label: "√3", observed: true, residual: 0.026 }], // 2.6%: past the 2.2% band, inside the 3% track
    };
    const { container } = render(<ResidualChart assigned={[s]} />);
    const pt = container.querySelector('[data-role="resid-point"]')!;
    expect(pt.getAttribute("data-outoftol")).toBe("true");
    expect(pt.querySelector("circle")!.getAttribute("fill")).toBe("none");       // hollow = out of tolerance
    expect(pt.querySelector('[data-role="resid-overflow"]')).toBeNull();          // still on-scale → not a chevron
  });

  it("positions the out-of-tolerance circle at its actual value, not clamped to the edge", () => {
    const mk = (residual: number): CombSeries => ({
      phase: "Pn3m", color: "oklch(0.570 0.150 58)",
      teeth: [{ q: 0.09, label: "√3", observed: true, residual }],
    });
    const a = render(<ResidualChart assigned={[mk(0.024)]} />);
    const b = render(<ResidualChart assigned={[mk(0.029)]} />);
    const cyA = Number(a.container.querySelector('[data-role="resid-point"] circle')!.getAttribute("cy"));
    const cyB = Number(b.container.querySelector('[data-role="resid-point"] circle')!.getAttribute("cy"));
    expect(cyB).toBeLessThan(cyA); // larger positive residual sits higher → distinct positions, not clamped to a shared edge
  });

  it("draws a clamped sign-flipped chevron and NO dot for a residual beyond the track", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    const offscale = container.querySelector('[data-role="resid-point"][data-overflow="true"]')!;
    expect(offscale.getAttribute("data-q")).toBe("0.1234"); // PN3M √6 residual +0.041 (> 3%)
    const arrow = offscale.querySelector('[data-role="resid-overflow"]')!;
    expect(arrow.getAttribute("data-dir")).toBe("up");      // positive → off the TOP
    expect(offscale.querySelector("circle")).toBeNull();    // off-scale → chevron only, no dot
  });

  it("points the off-scale chevron DOWN for a large negative residual", () => {
    const neg: CombSeries = {
      phase: "Pn3m", color: "oklch(0.570 0.150 58)",
      teeth: [{ q: 0.09, label: "√3", observed: true, residual: -0.05 }],
    };
    const { container } = render(<ResidualChart assigned={[neg]} />);
    const arrow = container.querySelector('[data-role="resid-overflow"]')!;
    expect(arrow.getAttribute("data-dir")).toBe("down");    // negative → off the BOTTOM
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
