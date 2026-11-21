import { describe, it, expect, vi } from "vitest";
import { render, fireEvent, screen } from "@testing-library/react";
import { COMPARE_COLORS, PLOT_HEIGHT, WorkspacePlot } from "../../src/workspace/ui/WorkspacePlot";
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
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [1, 0, 2], color: "var(--color-compare-1)", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    expect(d().match(/[ML]/g)).toHaveLength(2);
    rerender(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} view="linear" />);
    expect(d().match(/[ML]/g)).toHaveLength(3);
  });

  it("labels the q²I view and draws comparisons multiplied by q²", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [400, 100, 100], color: "var(--color-compare-1)", emphasis: false };
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
      .split(/[ML]/).filter(Boolean).map((p) => p.split(",").map(Number) as [number, number])
      .filter(([x]) => x >= 0 && x <= 640 - PLOT_MARGINS.left - PLOT_MARGINS.right);

  it("fits a stronger comparison into the window, after the view transform", () => {
    const flat = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const compare = { id: 3, name: "ref", q: flat.q, I: [5, 5, 5, 5, 5], color: "var(--color-compare-1)", emphasis: false };
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
    const short = { id: 4, name: "short", q: [0.02, 0.05, 0.08], I: [100, 100, 100], color: "var(--color-compare-1)", emphasis: false };
    render(<WorkspacePlot trace={flat} peaks={[]} combs={[]} compares={[short]} view="linear" xDomain={[0.32, 0.4]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
  });

  it("ignores a comparison that starts after the window", () => {
    const flat = { q: [0.05, 0.1, 0.15, 0.3, 0.4], I: [1, 1, 1, 1, 1], sigma: [0, 0, 0, 0, 0] };
    const late = { id: 5, name: "late", q: [0.3, 0.35, 0.4], I: [100, 100, 100], color: "var(--color-compare-1)", emphasis: false };
    render(<WorkspacePlot trace={flat} peaks={[]} combs={[]} compares={[late]} view="linear" xDomain={[0.05, 0.16]} />);
    expect(Math.max(...yTicks())).toBeLessThan(5);
  });

  it("falls back to the whole trace when the window holds no positive values", () => {
    const gap = { q: [0.05, 0.1, 0.12, 0.14, 0.2, 0.25], I: [40, 0, -3, Number.NaN, 0, 6], sigma: [0, 0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={gap} peaks={[]} combs={[]} view="linear" xDomain={[0.11, 0.15]} />);
    expect(Math.max(...yTicks())).toBeGreaterThanOrEqual(15);
  });

  it("ignores a strong peak just outside the window, so a zoom beside it is not flattened", () => {
    const spike = { q: [0.05, 0.1, 0.12, 0.2, 0.3, 0.4], I: [10, 100, 10, 10, 11, 10], sigma: [0, 0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={spike} peaks={[]} combs={[]} view="linear" xDomain={[0.15, 0.35]} />);
    expect(Math.max(...yTicks())).toBeLessThan(15);
  });

  it("keeps the line where it crosses the window edge inside the plot", () => {
    const ramp = { q: [0.05, 0.1, 0.2, 0.3, 0.4], I: [100, 50, 10, 10, 10], sigma: [0, 0, 0, 0, 0] };
    render(<WorkspacePlot trace={ramp} peaks={[]} combs={[]} view="linear" xDomain={[0.15, 0.35]} />);
    expect(Math.max(...yTicks())).toBeGreaterThanOrEqual(30);
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
    const glyph = (key: string): Element => screen.getByTestId("workspace-trace").querySelector(`[data-role="peak-glyph"] [data-peak-key="${key}"]`)!;
    const glyphYs = (key: string): number[] => glyph(key).getAttribute("points")!.split(" ").map((pt) => Number(pt.split(",")[1]));
    const { rerender } = render(<WorkspacePlot trace={steep} peaks={peaks} combs={[]} view="linear" />);
    const linear = lineYs();
    expect(linear[0]!).toBeLessThan(linear[linear.length - 1]!);
    const rowLinear = glyphYs("auto:1");
    rerender(<WorkspacePlot trace={steep} peaks={peaks} combs={[]} view="kratky" />);
    const kratky = lineYs();
    expect(kratky[0]!).toBeGreaterThan(kratky[kratky.length - 1]!);
    expect(glyphYs("auto:1")).toEqual(rowLinear);
    for (const key of ["auto:1", "auto:2", "manual:3"]) expect(Math.max(...glyphYs(key))).toBeLessThanOrEqual(0);
    expect(glyph("auto:1").getAttribute("fill")).toBe("var(--color-ink)");
    expect(glyph("auto:2").getAttribute("fill")).toBe("none");
    expect(glyph("auto:2").getAttribute("stroke")).toBe("var(--color-ink-faint)");
    expect(glyph("manual:3").getAttribute("points")!.split(" ")).toHaveLength(3);
    expect(glyph("manual:3").getAttribute("fill")).toBe("var(--color-accent)");
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

  const kpeaks = [
    { id: 1, exposure_id: 1, q: 0.1, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
    { id: 1, exposure_id: 1, q: 0.15, intensity: 8, prominence: null, sharpness: null, source: "manual" as const, excluded: true },
    { id: 2, exposure_id: 1, q: 0.3, intensity: 8, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
  ];

  it("lets the keyboard reach and toggle every peak in the window", () => {
    const onClickPeak = vi.fn();
    render(<WorkspacePlot trace={trace} peaks={kpeaks} combs={[]} onXDomain={() => {}} onClickPeak={onClickPeak} />);
    const auto = screen.getByRole("button", { name: "Auto peak at q = 0.1000 Å⁻¹" });
    const manual = screen.getByRole("button", { name: "Manual peak at q = 0.1500 Å⁻¹, off" });
    expect(auto.getAttribute("tabindex")).toBe("0");
    expect(screen.queryByRole("button", { name: /0\.3000/ })).toBeNull();
    fireEvent.keyDown(auto, { key: "Enter" });
    fireEvent.keyDown(manual, { key: " " });
    expect(onClickPeak.mock.calls).toEqual([[{ id: 1, source: "auto" }], [{ id: 1, source: "manual" }]]);
  });

  it("renders an auto and a manual peak that share an id without a duplicate key", () => {
    const error = vi.spyOn(console, "error").mockImplementation(() => {});
    render(<WorkspacePlot trace={trace} peaks={kpeaks} combs={[]} onXDomain={() => {}} onClickPeak={() => {}} />);
    expect(error.mock.calls.filter((c) => String(c[0]).includes("same key"))).toEqual([]);
    error.mockRestore();
  });

  it("a click reports the source of the peak it lands on when an auto and a manual peak share an id", () => {
    const onClickPeak = vi.fn();
    render(<WorkspacePlot trace={trace} peaks={kpeaks} combs={[]} onXDomain={() => {}} onClickPeak={onClickPeak} />);
    const plotWidth = 640 - PLOT_MARGINS.left - PLOT_MARGINS.right;
    const svg = screen.getByTestId("workspace-trace");
    fireEvent.click(svg, { clientX: PLOT_MARGINS.left + plotWidth * 0.5, clientY: 100 });
    fireEvent.click(svg, { clientX: PLOT_MARGINS.left + plotWidth * 0.25, clientY: 100 });
    expect(onClickPeak.mock.calls).toEqual([[{ id: 1, source: "manual" }], [{ id: 1, source: "auto" }]]);
  });

  it("does not hit a peak drawn outside the window, and adds peaks only inside the plot area", () => {
    const onClickPeak = vi.fn();
    const onAddPeak = vi.fn();
    const edge = [{ ...kpeaks[0]!, id: 7, q: 0.0995 }];
    render(<WorkspacePlot trace={trace} peaks={edge} combs={[]} xDomain={[0.1, 0.2]} onXDomain={() => {}} onClickPeak={onClickPeak} onAddPeak={onAddPeak} />);
    const svg = screen.getByTestId("workspace-trace");
    fireEvent.click(svg, { clientX: PLOT_MARGINS.left + 1, clientY: 100 });
    expect(onClickPeak).not.toHaveBeenCalled();
    expect(onAddPeak).toHaveBeenCalledTimes(1);
    fireEvent.click(svg, { clientX: PLOT_MARGINS.left + 200, clientY: 5 });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
    fireEvent.click(svg, { clientX: PLOT_MARGINS.left + 200, clientY: PLOT_HEIGHT - 5 });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
  });

  it("end-anchors an x label at the right edge only", () => {
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} xDomain={[0.1, 0.2]} />);
    const texts = Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-bottom"] text'));
    expect(texts.at(-1)!.getAttribute("text-anchor")).toBe("end");
    expect(texts[0]!.getAttribute("text-anchor")).toBeNull();
  });

  it("ignores pointer moves below the plot area", () => {
    if (typeof window.PointerEvent === "undefined") vi.stubGlobal("PointerEvent", MouseEvent);
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} />);
    const frame = screen.getByTestId("workspace-trace").parentElement!;
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + 100, clientY: PLOT_HEIGHT - 3 });
    expect(screen.getByTestId("workspace-trace").querySelector('[data-role="q-readout"]')).toBeNull();
  });

  it("gives each comparison its own colour with the mock-up's single dash", () => {
    const a = { id: 3, name: "a", q: trace.q, I: trace.I, color: COMPARE_COLORS[0], emphasis: true };
    const b = { id: 4, name: "b", q: trace.q, I: trace.I, color: COMPARE_COLORS[1], emphasis: false };
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[a, b]} />);
    const path = (id: number): Element => screen.getByTestId("workspace-compares").querySelector(`[data-compare="${id}"]`)!;
    expect(path(3).getAttribute("stroke")).not.toBe(path(4).getAttribute("stroke"));
    expect(path(3).getAttribute("stroke-dasharray")).toBe("5 4");
    expect(path(4).getAttribute("stroke-dasharray")).toBe("5 4");
    expect(path(3).getAttribute("stroke-width")).toBe("1.6");
  });

  it("draws the difference in the accent colour", () => {
    render(<SubtractionStrip q={[0.05, 0.1, 0.2]} d={[0, 10, 10]} label="minus buf" xDomain={null} xExtent={[0.05, 0.2]} peakQs={[]} />);
    const line = screen.getByTestId("workspace-subtraction").querySelector('[data-role="difference-line"]')!;
    expect(line.getAttribute("stroke")).toBe("var(--color-accent)");
    expect(line.getAttribute("stroke-width")).toBe("1.4");
  });
});
