import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CombChart, reflectionQExtent } from "../../src/print/comb/CombChart";
import type { CombSeries, CombRow } from "../../src/print/comb/combModel";

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

  it("dodges dense reflection labels apart with leader lines (FO-COMBLABEL-DODGE)", () => {
    // Three teeth packed into a tiny q-window → their centred labels would
    // overlap. They must dodge to >= the min gap and grow leader lines.
    const dense: CombSeries = {
      phase: "Im3m",
      color: "oklch(0.5 0.1 200)",
      teeth: [
        { q: 0.1, label: "√2", observed: true },
        { q: 0.101, label: "√4", observed: true },
        { q: 0.102, label: "√6", observed: true },
      ],
    };
    // Wide domain so the tight q-window maps to a tight pixel window (overlap).
    const { container } = render(<CombChart assigned={[dense]} leftover={[]} xDomain={[0.005, 0.3]} />);
    const xs = [...container.querySelectorAll('[data-role="tooth-mlabel"]')]
      .map((l) => Number(l.getAttribute("x")))
      .sort((a, b) => a - b);
    expect(xs.length).toBe(3);
    for (let i = 1; i < xs.length; i++) {
      expect(xs[i]! - xs[i - 1]!).toBeGreaterThanOrEqual(19.5); // ~LABEL_MIN_GAP
    }
    // shifted labels grow a leader back to their tooth
    expect(container.querySelectorAll('[data-role="tooth-leader"]').length).toBeGreaterThan(0);
  });

  it("does NOT draw leaders when the reflections are well separated", () => {
    const sparse: CombSeries = {
      phase: "Lamellar",
      color: "oklch(0.5 0.1 60)",
      teeth: [
        { q: 0.05, label: "1", observed: true },
        { q: 0.2, label: "2", observed: true },
      ],
    };
    const { container } = render(<CombChart assigned={[sparse]} leftover={[]} maxWidth={760} />);
    expect(container.querySelectorAll('[data-role="tooth-leader"]').length).toBe(0);
  });

  it("default-scrolls the q-pane to the reflections, not the left/beam edge (FO-COMBSCROLL-PEAKS)", () => {
    // A wide trace-linked domain → the comb overflows and would otherwise open at
    // scrollLeft 0 (the low-q beam dropoff). It must scroll right to the peaks.
    const { container } = render(
      <CombChart assigned={[PN3M]} leftover={[0.2]} xDomain={[0.005, 0.3]} />,
    );
    const pane = container.querySelector('[role="group"]') as HTMLDivElement;
    expect(pane).toBeTruthy();
    expect(pane.scrollLeft).toBeGreaterThan(0);
  });
});

describe("reflectionQExtent", () => {
  it("spans observed teeth + leftover rings, excluding predicted-absent carets", () => {
    const rows: CombRow[] = [
      { kind: "assigned", series: PN3M }, // observed 0.0712, 0.0872; absent 0.1126
      { kind: "leftover", qs: [0.2] },
    ];
    expect(reflectionQExtent(rows)).toEqual([0.0712, 0.2]);
  });

  it("falls back to all teeth when a row has only predicted-absent carets", () => {
    const absentOnly: CombSeries = {
      ...PN3M,
      teeth: [
        { q: 0.11, label: "a", observed: false },
        { q: 0.15, label: "b", observed: false },
      ],
    };
    expect(reflectionQExtent([{ kind: "assigned", series: absentOnly }])).toEqual([0.11, 0.15]);
  });

  it("returns undefined when there is nothing to scroll to", () => {
    expect(reflectionQExtent([])).toBeUndefined();
  });
});
