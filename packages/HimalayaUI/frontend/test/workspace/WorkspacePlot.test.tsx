import { describe, it, expect, vi } from "vitest";
import { render, fireEvent, screen } from "@testing-library/react";
import { PLOT_HEIGHT, WorkspacePlot } from "../../src/workspace/ui/WorkspacePlot";
import { PLOT_MARGINS } from "../../src/workspace/plotScale";
import { SubtractionStrip } from "../../src/workspace/ui/SubtractionStrip";

const trace = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [10, 50, 8, 6, 5], sigma: [0, 0, 0, 0, 0] };

describe("WorkspacePlot", () => {
  it("double-click does not reset the zoom", () => {
    const onXDomain = vi.fn();
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} xDomain={[0.1, 0.2]} onXDomain={onXDomain} />);
    fireEvent.doubleClick(screen.getByTestId("workspace-trace"));
    expect(onXDomain).not.toHaveBeenCalled();
  });

  it("skips non-positive comparison points on a log scale but draws them on linear", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [1, 0, 2], dash: "6 4", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    expect(d().match(/[ML]/g)).toHaveLength(2);
    rerender(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} view="linear" />);
    expect(d().match(/[ML]/g)).toHaveLength(3);
  });

  it("labels the q²I view and draws comparisons multiplied by q²", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [400, 100, 100], dash: "6 4", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} view="linear" />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    const ys = (): number[] => d().split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[1]));
    const linear = ys();
    expect(linear[0]!).toBeLessThan(linear[2]!);
    rerender(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} view="kratky" />);
    const kratky = ys();
    expect(kratky[0]!).toBeGreaterThan(kratky[2]!);
    expect(screen.getByTestId("workspace-trace").textContent).toContain("q²·I");
  });

  const linePoints = (): [number, number][] =>
    screen.getByTestId("workspace-trace").querySelector('[data-role="trace-line"] path:last-of-type')!.getAttribute("d")!
      .split(/[ML]/).filter(Boolean).map((p) => p.split(",").map(Number) as [number, number]);
  const yLabels = (): string[] =>
    Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-left"] text')).map((t) => t.textContent ?? "");
  const yTicks = (): number[] =>
    Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-left"] text'))
      .map((t) => Number(t.textContent))
      .filter(Number.isFinite);
  const plotHeight = (): number => PLOT_HEIGHT - PLOT_MARGINS.top - PLOT_MARGINS.bottom;

  it("fits the y axis to the visible q window in every view", () => {
    const rising = { q: [0.05, 0.1, 0.15, 0.2, 0.4], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const { rerender } = render(<WorkspacePlot trace={rising} peaks={[]} combs={[]} view="kratky" />);
    expect(Math.max(...yTicks())).toBeGreaterThan(0.1);
    rerender(<WorkspacePlot trace={rising} peaks={[]} combs={[]} view="kratky" xDomain={[0.1, 0.15]} />);
    const inWindow = linePoints().slice(1, 3).map(([, y]) => y);
    expect(Math.max(...inWindow) - Math.min(...inWindow)).toBeGreaterThan(plotHeight() * 0.25);
    const steep = { ...rising, I: [100, 50, 2, 1, 1] };
    rerender(<WorkspacePlot trace={steep} peaks={[]} combs={[]} view="linear" xDomain={[0.17, 0.4]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
    rerender(<WorkspacePlot trace={steep} peaks={[]} combs={[]} view="log" xDomain={[0.15, 0.4]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
    rerender(<WorkspacePlot trace={steep} peaks={[]} combs={[]} view="log" />);
    expect(yLabels()).toContain("10²");
  });

  it("keeps the tallest point below the top edge in the linear views", () => {
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} view="linear" />);
    const top = Math.min(...linePoints().map(([, y]) => y));
    expect(top).toBeGreaterThan(plotHeight() * 0.05);
  });

  const comparePoints = (): [number, number][] =>
    screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!
      .split(/[ML]/).filter(Boolean).map((p) => p.split(",").map(Number) as [number, number]);

  it("fits a stronger comparison into the window, after the view transform", () => {
    const flat = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const compare = { id: 3, name: "ref", q: flat.q, I: [5, 5, 5, 5, 5], dash: "6 4", emphasis: false };
    for (const view of ["linear", "kratky"] as const) {
      const { unmount } = render(<WorkspacePlot trace={flat} peaks={[]} combs={[]} compares={[compare]} view={view} xDomain={[0.1, 0.2]} />);
      const ys = comparePoints().map(([, y]) => y);
      expect(Math.min(...ys)).toBeGreaterThanOrEqual(0);
      expect(Math.min(...ys)).toBeLessThan(plotHeight() * 0.3);
      unmount();
    }
  });

  it("ignores a comparison that ends before the window", () => {
    const flat = { q: [0.05, 0.1, 0.3, 0.35, 0.4], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const short = { id: 4, name: "short", q: [0.02, 0.05, 0.08], I: [100, 100, 100], dash: "6 4", emphasis: false };
    render(<WorkspacePlot trace={flat} peaks={[]} combs={[]} compares={[short]} view="linear" xDomain={[0.32, 0.4]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
  });

  it("ignores a comparison that starts after the window", () => {
    const flat = { q: [0.05, 0.1, 0.15, 0.3, 0.4], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const late = { id: 5, name: "late", q: [0.3, 0.35, 0.4], I: [100, 100, 100], dash: "6 4", emphasis: false };
    render(<WorkspacePlot trace={flat} peaks={[]} combs={[]} compares={[late]} view="linear" xDomain={[0.05, 0.16]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
  });

  it("falls back to the whole trace when the window holds no positive values", () => {
    const gap = { q: [0.05, 0.1, 0.12, 0.14, 0.2, 0.25], I: [40, 0, -3, Number.NaN, 0, 6], sigma: [0, 0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={gap} peaks={[]} combs={[]} view="linear" xDomain={[0.11, 0.15]} />);
    expect(Math.max(...yTicks())).toBeGreaterThanOrEqual(15);
  });

  it("includes the samples just outside the window so the edge segments stay in range", () => {
    const ramp = { q: [0.05, 0.1, 0.2, 0.3, 0.4], I: [100, 50, 10, 10, 10], sigma: [0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={ramp} peaks={[]} combs={[]} view="linear" xDomain={[0.15, 0.35]} />);
    expect(Math.max(...yTicks())).toBeGreaterThanOrEqual(50);
  });

  it("includes the sample just after the window too", () => {
    const ramp = { q: [0.05, 0.1, 0.2, 0.3, 0.4], I: [10, 10, 10, 50, 100], sigma: [0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={ramp} peaks={[]} combs={[]} view="linear" xDomain={[0.08, 0.25]} />);
    expect(Math.max(...yTicks())).toBeGreaterThanOrEqual(50);
  });

  it("puts q on a log axis in the log-log view", () => {
    const xTicks = (): string =>
      Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-bottom"] text'))
        .map((t) => t.textContent)
        .join(" ");
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} view="log" />);
    const linearQ = xTicks();
    rerender(<WorkspacePlot trace={trace} peaks={[]} combs={[]} view="loglog" />);
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

  it("draws the trace multiplied by q² in the q²I view and keeps the peak markers in a row above the plot", () => {
    const steep = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [100, 40, 15, 8, 5], sigma: [0, 0, 0, 0, 0] };
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.05, intensity: 100, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
      { id: 2, exposure_id: 1, q: 0.25, intensity: 5, prominence: null, sharpness: null, source: "auto" as const, excluded: true },
      { id: 3, exposure_id: 1, q: 0.15, intensity: 15, prominence: null, sharpness: null, source: "manual" as const, excluded: false },
    ];
    const lineYs = (): number[] =>
      screen.getByTestId("workspace-trace").querySelector('[data-role="trace-line"] path')!.getAttribute("d")!
        .split(/[ML]/).filter(Boolean).map((p) => Number(p.split(",")[1]));
    const glyph = (id: number): Element => screen.getByTestId("workspace-trace").querySelector(`[data-role="peak-glyph"] [data-peak-id="${id}"]`)!;
    const glyphYs = (id: number): number[] => glyph(id).getAttribute("points")!.split(" ").map((pt) => Number(pt.split(",")[1]));
    const { rerender } = render(<WorkspacePlot trace={steep} peaks={peaks} combs={[]} view="linear" />);
    const linear = lineYs();
    expect(linear[0]!).toBeLessThan(linear[linear.length - 1]!);
    const rowLinear = glyphYs(1);
    rerender(<WorkspacePlot trace={steep} peaks={peaks} combs={[]} view="kratky" />);
    const kratky = lineYs();
    expect(kratky[0]!).toBeGreaterThan(kratky[kratky.length - 1]!);
    expect(glyphYs(1)).toEqual(rowLinear);
    for (const id of [1, 2, 3]) expect(Math.max(...glyphYs(id))).toBeLessThanOrEqual(0);
    expect(glyph(1).getAttribute("fill")).toBe("var(--color-ink)");
    expect(glyph(2).getAttribute("fill")).toBe("none");
    expect(glyph(3).getAttribute("points")!.split(" ")).toHaveLength(4);
  });

  it("shows q and d under the cursor and nothing once it leaves", () => {
    if (typeof window.PointerEvent === "undefined") vi.stubGlobal("PointerEvent", MouseEvent);
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} />);
    const frame = screen.getByTestId("workspace-trace").parentElement!;
    const plotWidth = 640 - PLOT_MARGINS.left - PLOT_MARGINS.right;
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + plotWidth * 0.25, clientY: 100 });
    expect(screen.getByTestId("workspace-trace").querySelector('[data-role="q-readout"]')!.textContent).toBe("q 0.1000 Å⁻¹ · d 62.8 Å");
    fireEvent.pointerLeave(frame);
    expect(screen.getByTestId("workspace-trace").querySelector('[data-role="q-readout"]')).toBeNull();
  });
});
