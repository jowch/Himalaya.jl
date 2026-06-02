import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { WaterfallChart } from "../../src/print/waterfall/WaterfallChart";
import type { WaterfallRow } from "../../src/print/waterfall/waterfallModel";

const ROWS: WaterfallRow[] = [
  {
    key: "a", label: "1:0", phase: "Ia3d", state: "indexed",
    trace: { q: [0.02, 0.05, 0.1], I: [100, 60, 12], sigma: [1, 1, 1] },
    anchors: [{ id: 1, q: 0.05, intensity: 60, phase: "Ia3d" }],
    bandHeight: 1, yOffset: 0,
  },
  {
    key: "b", label: "1:1", phase: "Im3m", state: "indexed",
    trace: { q: [0.02, 0.06, 0.12], I: [90, 55, 10], sigma: [1, 1, 1] },
    anchors: [{ id: 2, q: 0.06, intensity: 55, phase: "Im3m" }],
    bandHeight: 1, yOffset: 0,
  },
];

describe("WaterfallChart", () => {
  it("renders one row group per member", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    expect(container.querySelectorAll('[data-role="wf-row"]').length).toBe(2);
  });

  it("stacks display-order-0 at the BOTTOM (largest top offset)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} maxWidth={800} />);
    const groups = Array.from(container.querySelectorAll('[data-role="wf-row"]'));
    const topA = Number((groups.find((g) => g.getAttribute("data-key") === "a") as HTMLElement).style.top.replace("px", ""));
    const topB = Number((groups.find((g) => g.getAttribute("data-key") === "b") as HTMLElement).style.top.replace("px", ""));
    expect(topA).toBeGreaterThan(topB);
  });

  it("renders a phase-coloured trace line per row (not ink-soft)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    const lines = container.querySelectorAll('[data-role="trace-line"] path');
    expect(lines.length).toBeGreaterThanOrEqual(2);
    const strokes = Array.from(lines).map((p) => p.getAttribute("stroke"));
    expect(strokes.some((s) => s?.includes("oklch"))).toBe(true);
  });

  it("renders a bead per anchor (the peaks layer)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    expect(container.querySelectorAll('[data-role="peak-glyph"]').length).toBeGreaterThanOrEqual(2);
  });

  it("renders exactly one shared bottom axis", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    expect(container.querySelectorAll('[data-role="axis-bottom"]').length).toBe(1);
  });

  it("renders the sample label per row", () => {
    const { getByText } = render(<WaterfallChart rows={ROWS} />);
    expect(getByText("1:0")).toBeTruthy();
    expect(getByText("1:1")).toBeTruthy();
  });

  it("marks the hovered row hot and dims the others; fires onHoverRow", () => {
    const spy = vi.fn();
    const { container } = render(<WaterfallChart rows={ROWS} onHoverRow={spy} />);
    const rowA = container.querySelector('[data-role="wf-row"][data-key="a"]')!;
    fireEvent.mouseEnter(rowA);
    expect(spy).toHaveBeenCalledWith("a");
    const rowB = container.querySelector('[data-role="wf-row"][data-key="b"]')!;
    expect(rowA.getAttribute("data-hot")).toBe("true");
    expect(rowB.getAttribute("data-dim")).toBe("true");
    fireEvent.mouseLeave(rowA);
    expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("respects a controlled hoveredKey", () => {
    const { container } = render(<WaterfallChart rows={ROWS} hoveredKey="b" />);
    expect(container.querySelector('[data-role="wf-row"][data-key="b"]')!.getAttribute("data-hot")).toBe("true");
    expect(container.querySelector('[data-role="wf-row"][data-key="a"]')!.getAttribute("data-dim")).toBe("true");
  });

  it("defaults the q-axis to log and exposes it via data-xtype", () => {
    const { getByTestId } = render(<WaterfallChart rows={ROWS} />);
    expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("log");
  });

  it("renders a log/linear scale toggle and flips internal scale on click", () => {
    const { getByTestId } = render(<WaterfallChart rows={ROWS} />);
    expect(getByTestId("wf-scale")).toBeTruthy();
    const linear = getByTestId("wf-scale").querySelector('[data-value="linear"]') as HTMLElement;
    fireEvent.click(linear);
    expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("linear");
  });

  it("fires onXTypeChange and respects a controlled xType", () => {
    const onXTypeChange = vi.fn();
    const { getByTestId, rerender } = render(
      <WaterfallChart rows={ROWS} xType="log" onXTypeChange={onXTypeChange} />,
    );
    const linear = getByTestId("wf-scale").querySelector('[data-value="linear"]') as HTMLElement;
    fireEvent.click(linear);
    expect(onXTypeChange).toHaveBeenCalledWith("linear");
    expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("log");
    rerender(<WaterfallChart rows={ROWS} xType="linear" onXTypeChange={onXTypeChange} />);
    expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("linear");
  });

  it("renders no q-guide when hoveredQ is unset", () => {
    const { queryByTestId } = render(<WaterfallChart rows={ROWS} />);
    expect(queryByTestId("wf-qguide")).toBeNull();
    expect(queryByTestId("wf-qreadout")).toBeNull();
  });

  it("renders the q-guide + readout at a controlled hoveredQ", () => {
    const q = ROWS[0]!.anchors[0]!.q;
    const { getByTestId } = render(<WaterfallChart rows={ROWS} hoveredQ={q} />);
    expect(getByTestId("wf-qguide").getAttribute("data-q")).toBe(String(q));
    expect(getByTestId("wf-qreadout").textContent).toBe(q.toFixed(3));
  });

  it("clears the cursor (fires onHoverQ undefined) on pointer leave", () => {
    const onHoverQ = vi.fn();
    const { getByTestId } = render(<WaterfallChart rows={ROWS} onHoverQ={onHoverQ} />);
    const stack = getByTestId("wf-stack");
    fireEvent.pointerLeave(stack);
    expect(onHoverQ).toHaveBeenLastCalledWith(undefined);
  });
});
