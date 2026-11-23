import { describe, it, expect, vi } from "vitest";
import { act, render, fireEvent, screen } from "@testing-library/react";
import { COMPARE_COLORS, PLOT_HEIGHT, WorkspacePlot } from "../../src/ui/WorkspacePlot";
import { PLOT_MARGINS } from "../../src/workspace/plotScale";

const trace = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [10, 50, 8, 6, 5], sigma: [0, 0, 0, 0, 0] };

describe("WorkspacePlot", () => {
  it("says so when a frame's trace has no points, but stays quiet while it loads", () => {
    const { rerender } = render(<WorkspacePlot trace={undefined} peaks={[]} combs={[]} />);
    expect(screen.getByTestId("workspace-trace-empty").textContent).toBe("");
    rerender(<WorkspacePlot trace={{ q: [], I: [], sigma: [] } as never} peaks={[]} combs={[]} />);
    expect(screen.getByTestId("workspace-trace-empty")).toHaveTextContent("This image's trace has no points.");
  });

  it("puts the q label in the tick row, below the data, without a tick label under it", () => {
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} />);
    const svg = screen.getByTestId("workspace-trace").closest("div")!.querySelector("svg")!;
    const label = [...svg.querySelectorAll("text")].find((t) => t.textContent === "q (Å⁻¹)")!;
    const plotHeight = PLOT_HEIGHT - PLOT_MARGINS.top - PLOT_MARGINS.bottom;
    expect(Number(label.getAttribute("y"))).toBeGreaterThan(plotHeight);
    const labelX = Number(label.getAttribute("x"));
    const ticks = [...svg.querySelectorAll('[data-role="axis-bottom"] text')].filter((t) => t !== label);
    expect(ticks.length).toBeGreaterThan(0);
    for (const t of ticks) expect(Number(t.getAttribute("x"))).toBeLessThan(labelX - 56);
  });

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

  const diffPoints = (): [number, number][] =>
    screen.getByTestId("workspace-trace").querySelector('[data-role="difference-line"] path')!.getAttribute("d")!
      .split(/[ML]/).filter(Boolean).map((p) => p.split(",").map(Number) as [number, number]);

  it("draws the difference in place of the trace, with a zero line and no comparisons", () => {
    const compare = { id: 3, name: "buf", q: trace.q, I: trace.I, color: COMPARE_COLORS[0], emphasis: true };
    const difference = { q: [0.05, 0.1, 0.15, 0.2, 0.25], d: [4, 40, -2, 1, -1] };
    render(<WorkspacePlot trace={trace} peaks={[]} combs={[]} compares={[compare]} difference={difference} view="linear" />);
    const svg = screen.getByTestId("workspace-trace");
    expect(svg.querySelector('[data-role="trace-line"]')).toBeNull();
    expect(svg.querySelector('[data-role="difference-line"] path')!.getAttribute("stroke")).toBe("var(--color-accent)");
    expect(screen.getByTestId("workspace-compares").querySelector("path")).toBeNull();
    const zero = Number(svg.querySelector('[data-role="zero-line"]')!.getAttribute("y1"));
    const pts = diffPoints();
    expect(pts[2]![1]).toBeGreaterThan(zero);
    expect(pts[3]![1]).toBeLessThan(zero);
    expect(yTicks()).toContain(0);
  });

  it("puts the difference on a pseudolog axis in the log views and a linear one otherwise", () => {
    const q = Array.from({ length: 40 }, (_, i) => 0.05 + i * 0.005);
    const d = q.map((v, i) => (Math.abs(v - 0.12) < 0.004 ? 1000 : i % 2 ? 1 : -1));
    const props = { trace, peaks: [], combs: [], difference: { q, d } };
    const { rerender } = render(<WorkspacePlot {...props} view="linear" />);
    expect(yLabels().some((t) => t.startsWith("−") || t.startsWith("-"))).toBe(true);
    const linear = yTicks();
    rerender(<WorkspacePlot {...props} view="log" />);
    const pseudo = yTicks();
    expect(pseudo).toContain(0);
    expect(pseudo).not.toEqual(linear);
    expect(pseudo.filter((t) => t > 0).every((t) => Math.log10(t) % 1 === 0)).toBe(true);
    rerender(<WorkspacePlot {...props} view="kratky" />);
    expect(yTicks()).not.toEqual(pseudo);
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

  it("draws a dashed crosshair at the readout's q, under the peak markers, in every view and zoom", () => {
    if (typeof window.PointerEvent === "undefined") vi.stubGlobal("PointerEvent", MouseEvent);
    const peaks = [{ id: 1, exposure_id: 1, q: 0.1, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false }];
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} />);
    const svg = (): Element => screen.getByTestId("workspace-trace");
    const crosshair = (): Element | null => svg().querySelector('[data-role="crosshair"]');
    expect(crosshair()).toBeNull();
    const frame = svg().parentElement!;
    const plotWidth = 640 - PLOT_MARGINS.left - PLOT_MARGINS.right;
    const at = plotWidth * 0.25;
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + at, clientY: 100 });
    const line = crosshair()!;
    expect(Number(line.getAttribute("x1"))).toBeCloseTo(at, 6);
    expect(line.getAttribute("x2")).toBe(line.getAttribute("x1"));
    expect(Number(line.getAttribute("y1"))).toBe(0);
    expect(Number(line.getAttribute("y2"))).toBe(PLOT_HEIGHT - PLOT_MARGINS.top - PLOT_MARGINS.bottom);
    expect(line.getAttribute("stroke-dasharray")).toBeTruthy();
    expect(line.closest('[pointer-events="none"]') ?? (line.getAttribute("pointer-events") === "none" ? line : null)).not.toBeNull();
    const glyphs = svg().querySelector('[data-role="peak-glyph"]')!;
    expect(line.compareDocumentPosition(glyphs) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    for (const view of ["linear", "loglog", "kratky"] as const) {
      rerender(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} view={view} />);
      fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + at, clientY: 100 });
      const readout = Number(svg().querySelector('[data-role="q-readout"]')!.textContent!.match(/q ([\d.]+)/)![1]);
      expect(readout).toBeGreaterThan(0);
      expect(Number(crosshair()!.getAttribute("x1"))).toBeCloseTo(at, 6);
    }
    rerender(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} xDomain={[0.1, 0.2]} />);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + at, clientY: 100 });
    expect(svg().querySelector('[data-role="q-readout"]')!.textContent).toMatch(/^q 0\.12/);
    expect(Number(crosshair()!.getAttribute("x1"))).toBeCloseTo(at, 6);
    fireEvent.pointerLeave(frame);
    expect(crosshair()).toBeNull();
  });

  it("draws a faded comb lighter than a hovered one, at the same q", () => {
    const combs = [
      { id: 1, qs: [0.1], color: "red" },
      { id: 2, qs: [0.1], color: "blue", faded: true },
    ];
    render(<WorkspacePlot trace={trace} peaks={[]} combs={combs} />);
    const svg = screen.getByTestId("workspace-trace");
    const full = svg.querySelector('[data-comb="1"]')!;
    const faded = svg.querySelector('[data-comb="2"]')!;
    expect(faded.getAttribute("x1")).toBe(full.getAttribute("x1"));
    expect(full.getAttribute("stroke-opacity")).toBeNull();
    expect(Number(faded.getAttribute("stroke-opacity"))).toBeGreaterThan(0);
    expect(Number(faded.getAttribute("stroke-opacity"))).toBeLessThan(0.5);
  });

  it("outlines an unlit manual marker in its own colour so it reads as large as an auto marker on a dark plot", () => {
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.1, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
      { id: 2, exposure_id: 1, q: 0.15, intensity: 50, prominence: null, sharpness: null, source: "manual" as const, excluded: false },
      { id: 3, exposure_id: 1, q: 0.2, intensity: 50, prominence: null, sharpness: null, source: "manual" as const, excluded: true },
    ];
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} />);
    const glyph = (key: string): Element => screen.getByTestId("workspace-trace").querySelector(`[data-peak-key="${key}"]`)!;
    expect(glyph("manual:2").getAttribute("stroke")).toBe("var(--color-accent)");
    expect(glyph("manual:2").classList.contains("peak-manual")).toBe(true);
    expect(glyph("auto:1").classList.contains("peak-manual")).toBe(false);
    expect(glyph("manual:3").classList.contains("peak-manual")).toBe(false);
    rerender(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} highlight={new Map([["manual:2", "red"]])} />);
    expect(glyph("manual:2").classList.contains("peak-manual")).toBe(false);
    expect(glyph("manual:2").getAttribute("stroke")).toBe("red");
  });

  it("snaps the crosshair and readout to a peak's q while the pointer is on its marker, and not in the plot below it", () => {
    if (typeof window.PointerEvent === "undefined") vi.stubGlobal("PointerEvent", MouseEvent);
    const peaks = [{ id: 1, exposure_id: 1, q: 0.1234, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false }];
    render(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} />);
    const svg = screen.getByTestId("workspace-trace");
    const frame = svg.parentElement!;
    const tip = Number(svg.querySelector('[data-peak-key="auto:1"]')!.getAttribute("points")!.split(" ")[2]!.split(",")[0]);
    const readout = (): string => svg.querySelector('[data-role="q-readout"]')!.textContent ?? "";
    const line = (): number => Number(svg.querySelector('[data-role="crosshair"]')!.getAttribute("x1"));
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + tip + 3, clientY: PLOT_MARGINS.top - 4 });
    expect(readout()).toMatch(/^q 0\.1234 /);
    expect(line()).toBeCloseTo(tip, 1);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + tip + 3, clientY: 100 });
    expect(readout()).not.toMatch(/^q 0\.1234 /);
    expect(line()).toBeCloseTo(tip + 3, 1);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + tip + 3, clientY: 2 });
    expect(readout()).not.toMatch(/^q 0\.1234 /);
    expect(line()).toBeCloseTo(tip + 3, 1);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + tip + 40, clientY: PLOT_MARGINS.top - 4 });
    expect(line()).toBeCloseTo(tip + 40, 1);
  });

  it("shows a crosshair cursor exactly while the crosshair is drawn, and keeps the pointer on peak markers", () => {
    if (typeof window.PointerEvent === "undefined") vi.stubGlobal("PointerEvent", MouseEvent);
    const peaks = [{ id: 1, exposure_id: 1, q: 0.1, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false }];
    render(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} onClickPeak={() => {}} onAddPeak={() => {}} onXDomain={() => {}} />);
    const svg = screen.getByTestId("workspace-trace");
    const frame = svg.parentElement!;
    const plotWidth = 640 - PLOT_MARGINS.left - PLOT_MARGINS.right;
    expect(svg.classList.contains("cursor-crosshair")).toBe(false);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + plotWidth * 0.25, clientY: 100 });
    expect(svg.querySelector('[data-role="crosshair"]')).not.toBeNull();
    expect(svg.classList.contains("cursor-crosshair")).toBe(true);
    expect(svg.querySelector("[data-peak-key]")!.classList.contains("cursor-pointer")).toBe(true);
    fireEvent.pointerMove(frame, { clientX: 2, clientY: 100 });
    expect(svg.querySelector('[data-role="crosshair"]')).toBeNull();
    expect(svg.classList.contains("cursor-crosshair")).toBe(false);
    fireEvent.pointerMove(frame, { clientX: PLOT_MARGINS.left + plotWidth * 0.5, clientY: 100 });
    fireEvent.pointerLeave(frame);
    expect(svg.classList.contains("cursor-crosshair")).toBe(false);
  });

  it("colours the peaks a hovered call claims and dims the rest, without moving any marker", () => {
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.1, intensity: 50, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
      { id: 2, exposure_id: 1, q: 0.15, intensity: 8, prominence: null, sharpness: null, source: "auto" as const, excluded: false },
      { id: 3, exposure_id: 1, q: 0.2, intensity: 6, prominence: null, sharpness: null, source: "manual" as const, excluded: false },
    ];
    const glyph = (key: string): Element => screen.getByTestId("workspace-trace").querySelector(`[data-peak-key="${key}"]`)!;
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} />);
    const before = ["auto:1", "auto:2", "manual:3"].map((k) => glyph(k).getAttribute("points"));
    for (const k of ["auto:1", "auto:2", "manual:3"]) expect(glyph(k).getAttribute("opacity")).toBeNull();
    rerender(<WorkspacePlot trace={trace} peaks={peaks} combs={[]} highlight={new Map([["auto:1", "var(--color-phase-pn3m)"], ["manual:3", "var(--color-phase-im3m)"]])} />);
    expect(glyph("auto:1").getAttribute("fill")).toBe("var(--color-phase-pn3m)");
    expect(glyph("manual:3").getAttribute("fill")).toBe("var(--color-phase-im3m)");
    expect(glyph("auto:1").getAttribute("opacity")).toBeNull();
    expect(Number(glyph("auto:2").getAttribute("opacity"))).toBeLessThan(1);
    expect(glyph("auto:2").getAttribute("fill")).toBe("var(--color-ink)");
    expect(["auto:1", "auto:2", "manual:3"].map((k) => glyph(k).getAttribute("points"))).toEqual(before);
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

  it("a pointer press on a peak does not take keyboard focus, so Enter still reaches the workspace", () => {
    render(<WorkspacePlot trace={trace} peaks={kpeaks} combs={[]} onXDomain={() => {}} onClickPeak={() => {}} />);
    const auto = screen.getByRole("button", { name: "Auto peak at q = 0.1000 Å⁻¹" });
    expect(fireEvent.mouseDown(auto)).toBe(false);
  });

  it("rings the keyboard-focused peak without recolouring it", () => {
    render(<WorkspacePlot trace={trace} peaks={kpeaks} combs={[]} onXDomain={() => {}} onClickPeak={() => {}} />);
    const svg = screen.getByTestId("workspace-trace");
    const ring = (): Element | null => svg.querySelector('[data-role="peak-focus"]');
    const auto = screen.getByRole("button", { name: "Auto peak at q = 0.1000 Å⁻¹" });
    expect(ring()).toBeNull();
    act(() => (auto as unknown as HTMLElement).focus());
    expect(ring()).not.toBeNull();
    expect(ring()!.getAttribute("fill")).toBe("none");
    expect(["var(--color-accent)", "var(--color-ink)"]).not.toContain(ring()!.getAttribute("stroke"));
    expect(auto.getAttribute("fill")).toBe("var(--color-ink)");
    expect(auto.getAttribute("stroke")).toBe("none");
    act(() => (auto as unknown as HTMLElement).blur());
    expect(ring()).toBeNull();
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

});
