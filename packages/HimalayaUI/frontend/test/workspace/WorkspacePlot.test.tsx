import { describe, it, expect, vi } from "vitest";
import { render, fireEvent, screen } from "@testing-library/react";
import { WorkspacePlot } from "../../src/workspace/ui/WorkspacePlot";
import { SubtractionStrip } from "../../src/workspace/ui/SubtractionStrip";

const trace = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [10, 50, 8, 6, 5], sigma: [0, 0, 0, 0, 0] };

describe("WorkspacePlot", () => {
  it("double-click does not reset the zoom", () => {
    const onXDomain = vi.fn();
    render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} xDomain={[0.1, 0.2]} onXDomain={onXDomain} />);
    fireEvent.doubleClick(screen.getByTestId("workspace-trace"));
    expect(onXDomain).not.toHaveBeenCalled();
  });

  it("skips non-positive comparison points on a log scale but draws them on linear", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [1, 0, 2], dash: "6 4", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    expect(d().match(/[ML]/g)).toHaveLength(2);
    rerender(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} view="linear" />);
    expect(d().match(/[ML]/g)).toHaveLength(3);
  });

  it("labels the q²I view and draws comparisons multiplied by q²", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [400, 100, 100], dash: "6 4", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} view="linear" />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    const ys = (): number[] => d().split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[1]));
    const linear = ys();
    expect(linear[0]!).toBeLessThan(linear[2]!);
    rerender(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} view="kratky" />);
    const kratky = ys();
    expect(kratky[0]!).toBeGreaterThan(kratky[2]!);
    expect(screen.getByTestId("workspace-trace").textContent).toContain("q²·I");
  });

  it("fits the q²I and linear y axes to the visible q window, but not the log axes", () => {
    const rising = { q: [0.05, 0.1, 0.15, 0.2, 0.4], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const yTop = (): number =>
      Math.max(
        ...Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-left"] text'))
          .map((t) => Number(t.textContent))
          .filter(Number.isFinite),
      );
    const { rerender } = render(<WorkspacePlot trace={rising} peaks={[]} phase={null} combs={[]} view="kratky" />);
    expect(yTop()).toBeGreaterThan(0.1);
    rerender(<WorkspacePlot trace={rising} peaks={[]} phase={null} combs={[]} view="kratky" xDomain={[0.05, 0.15]} />);
    expect(yTop()).toBeLessThan(0.03);
    const steep = { ...rising, I: [100, 50, 2, 1, 1] };
    rerender(<WorkspacePlot trace={steep} peaks={[]} phase={null} combs={[]} view="linear" xDomain={[0.15, 0.4]} />);
    expect(yTop()).toBeLessThan(5);
    rerender(<WorkspacePlot trace={steep} peaks={[]} phase={null} combs={[]} view="log" xDomain={[0.15, 0.4]} />);
    expect(yTop()).toBeGreaterThanOrEqual(100);
  });

  it("puts q on a log axis in the log-log view", () => {
    const xTicks = (): string =>
      Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-bottom"] text'))
        .map((t) => t.textContent)
        .join(" ");
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} view="log" />);
    const linearQ = xTicks();
    rerender(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} view="loglog" />);
    expect(xTicks()).not.toBe(linearQ);
    expect(xTicks().split(" ")).toContain("0.10");
  });

  it("keeps the difference strip on the plot's q axis and intensity transform", () => {
    const props = { q: [0.05, 0.1, 0.2], d: [0, 10, 10], label: "minus buf", xDomain: null, xExtent: [0.05, 0.2] as [number, number], peakQs: [] };
    const xs = (): number[] =>
      screen.getByTestId("workspace-subtraction").querySelector("[data-role='difference-line']")!.getAttribute("d")!
        .split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[0]));
    const ys = (): number[] =>
      screen.getByTestId("workspace-subtraction").querySelector("[data-role='difference-line']")!.getAttribute("d")!
        .split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[1]));
    const { rerender } = render(<SubtractionStrip {...props} view="log" />);
    const linearX = xs();
    expect(ys()[1]).toBeCloseTo(ys()[2]!);
    rerender(<SubtractionStrip {...props} view="loglog" />);
    expect(xs()[1]! - xs()[0]!).toBeGreaterThan(linearX[1]! - linearX[0]!);
    rerender(<SubtractionStrip {...props} view="kratky" />);
    expect(ys()[2]!).toBeLessThan(ys()[1]!);
  });

  it("draws the trace and its peak markers multiplied by q² in the q²I view", () => {
    const steep = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [100, 40, 15, 8, 5], sigma: [0, 0, 0, 0, 0] };
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.05, intensity: 100, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
      { id: 2, exposure_id: 1, q: 0.25, intensity: 5, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
    ];
    const lineYs = (): number[] =>
      screen.getByTestId("workspace-trace").querySelector('[data-role="trace-line"] path:last-of-type')!.getAttribute("d")!
        .split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[1]));
    const markerY = (id: number): number => {
      const glyph = screen.getByTestId("workspace-trace").querySelector(`[data-role="peak-glyph"] [data-peak-id="${id}"]`)!;
      return Math.max(...glyph.getAttribute("points")!.split(" ").map((pt) => Number(pt.split(",")[1])));
    };
    const { rerender } = render(<WorkspacePlot trace={steep} peaks={peaks} phase={null} combs={[]} view="linear" />);
    const linear = lineYs();
    expect(linear[0]!).toBeLessThan(linear[linear.length - 1]!);
    expect(markerY(1)).toBeLessThan(markerY(2));
    rerender(<WorkspacePlot trace={steep} peaks={peaks} phase={null} combs={[]} view="kratky" />);
    const kratky = lineYs();
    expect(kratky[0]!).toBeGreaterThan(kratky[kratky.length - 1]!);
    expect(markerY(1)).toBeGreaterThan(markerY(2));
  });
});
