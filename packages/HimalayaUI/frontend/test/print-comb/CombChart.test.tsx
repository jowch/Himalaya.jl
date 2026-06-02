import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CombChart } from "../../src/print/comb/CombChart";
import type { CombSeries } from "../../src/print/comb/combModel";

const PN3M: CombSeries = {
  phase: "Pn3m",
  color: "oklch(0.570 0.150 58)",
  latticeLabel: "a = 197 Å",
  teeth: [
    { q: 0.0712, label: "√2", observed: true },
    { q: 0.0872, label: "√3", observed: true },
    { q: 0.1126, label: "√5", observed: false }, // predicted-absent
  ],
};

describe("CombChart", () => {
  it("renders one comb-row per assigned series plus a leftover row", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[0.2]} />);
    expect(container.querySelectorAll('[data-role="comb-row"]').length).toBe(2);
    expect(container.querySelector('[data-row-kind="leftover"]')).toBeTruthy();
  });

  it("draws a solid tooth (stem + cap + label) for observed reflections", () => {
    const { container, getByText } = render(<CombChart assigned={[PN3M]} leftover={[]} />);
    expect(container.querySelectorAll('[data-role="tooth"]').length).toBe(2);
    expect(container.querySelector('[data-role="tooth-cap"]')).toBeTruthy();
    expect(getByText("√2")).toBeTruthy();
    expect(container.querySelector('[data-role="tooth-stem"]')!.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
  });

  it("draws a dashed caret (not a tooth) for a predicted-absent reflection", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} />);
    expect(container.querySelectorAll('[data-role="caret"]').length).toBe(1);
    const stem = container.querySelector('[data-role="caret-stem"]')!;
    expect(stem.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("draws a ring per leftover peak", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[0.2, 0.25]} />);
    expect(container.querySelectorAll('[data-role="leftover-ring"]').length).toBe(2);
  });

  it("renders a dashed preview row when a hovered candidate is supplied", () => {
    const { container } = render(<CombChart assigned={[PN3M]} hovered={PN3M} leftover={[]} />);
    expect(container.querySelector('[data-role="row-baseline"][data-preview="true"]')).toBeTruthy();
  });

  it("a hovered tooth goes hot, keeping its OWN colour (no accent recolor) + a ring", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} hoveredQ={0.0712} />);
    const hot = container.querySelector('[data-role="tooth"][data-hot="true"]')!;
    expect(hot.getAttribute("data-q")).toBe("0.0712");
    const stem = hot.querySelector('[data-role="tooth-stem"]')!;
    expect(stem.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
    expect(stem.getAttribute("stroke")).not.toBe("var(--color-accent)");
    expect(hot.querySelector('[data-role="tooth-ring"]')).toBeTruthy();
  });

  it("fires onHoverQ(q) on tooth enter and onHoverQ(undefined) on leave", () => {
    const spy = vi.fn();
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} onHoverQ={spy} />);
    const tooth = container.querySelector('[data-role="tooth"]')!;
    fireEvent.mouseEnter(tooth);
    expect(spy).toHaveBeenCalledWith(0.0712);
    fireEvent.mouseLeave(tooth);
    expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("enforces the min plot width at a narrow maxWidth (pane scrolls)", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} maxWidth={320} />);
    expect(Number(container.querySelector("svg")!.getAttribute("data-plot-w"))).toBeGreaterThanOrEqual(320);
  });
});
