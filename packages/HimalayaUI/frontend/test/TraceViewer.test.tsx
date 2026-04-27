import { describe, it, expect, vi } from "vitest";
import { render } from "@testing-library/react";
import { TraceViewer, findNearestPeak } from "../src/components/TraceViewer";
import type { IndexEntry } from "../src/api";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    // Return a plot-like element with a stub `.scale` so the component's
    // overlay code path doesn't crash.
    const el = document.createElement("div");
    el.setAttribute("data-testid", "plot-svg");
    (el as unknown as { scale: (n: string) => unknown }).scale = () => ({
      apply: (q: number) => q * 1000,
      invert: (px: number) => px / 1000,
    });
    return el;
  }),
  areaY: vi.fn(() => ({ _kind: "areaY" })),
  line:  vi.fn(() => ({ _kind: "line" })),
  dot:   vi.fn(() => ({ _kind: "dot" })),
}));

describe("findNearestPeak", () => {
  const peaks = [
    { id: 1, exposure_id: 1, q: 0.10, intensity: null, prominence: null,
      sharpness: null, source: "auto" as const, excluded: false },
    { id: 2, exposure_id: 1, q: 0.30, intensity: null, prominence: null,
      sharpness: null, source: "auto" as const, excluded: false },
  ];

  it("returns the closest peak when within tolerance", () => {
    expect(findNearestPeak(peaks, 0.105, 0.01)?.id).toBe(1);
  });

  it("returns null when outside tolerance", () => {
    expect(findNearestPeak(peaks, 0.20, 0.01)).toBeNull();
  });

  it("returns null for empty list", () => {
    expect(findNearestPeak([], 0.1, 1)).toBeNull();
  });
});

const defaultProps = {
  onAddPeak: () => {},
  onRemovePeak: () => {},
  onTogglePeakExclusion: () => {},
  xDomain: null as [number, number] | null,
  onXDomain: () => {},
};

describe("<TraceViewer>", () => {
  it("renders a container with testid 'trace-viewer'", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        activeGroupIndices={[]}
        hoveredIndex={undefined}
        {...defaultProps}
      />,
    );
    expect(container.querySelector('[data-testid="trace-viewer"]')).not.toBeNull();
  });

  it("invokes Plot.plot with marks including line and areaY for trace", async () => {
    const Plot = await import("@observablehq/plot");
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
    render(
      <TraceViewer trace={trace} peaks={[]} activeGroupIndices={[]} hoveredIndex={undefined}
        {...defaultProps} />,
    );
    expect(Plot.plot).toHaveBeenCalled();
    expect(Plot.areaY).toHaveBeenCalled();
    expect(Plot.line).toHaveBeenCalled();
  });

  it("renders one overlay triangle per peak (auto + manual together)", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.1, intensity: null, prominence: null,
        sharpness: null, source: "auto"   as const, excluded: false },
      { id: 2, exposure_id: 1, q: 0.2, intensity: null, prominence: null,
        sharpness: null, source: "manual" as const, excluded: false },
    ];
    const { container } = render(
      <TraceViewer trace={trace} peaks={peaks} activeGroupIndices={[]} hoveredIndex={undefined}
        {...defaultProps} />,
    );
    const triangles = container.querySelectorAll('[data-role="peak-root"] polygon');
    expect(triangles.length).toBe(2);
  });

  it("excluded auto peaks render at reduced opacity", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const peaks = [
      { id: 1, exposure_id: 1, q: 0.1, intensity: null, prominence: null,
        sharpness: null, source: "auto"   as const, excluded: true },
      { id: 2, exposure_id: 1, q: 0.2, intensity: null, prominence: null,
        sharpness: null, source: "auto"   as const, excluded: false },
    ];
    const { container } = render(
      <TraceViewer trace={trace} peaks={peaks} activeGroupIndices={[]} hoveredIndex={undefined}
        {...defaultProps} />,
    );
    const triangles = container.querySelectorAll('[data-role="peak-root"] polygon');
    const opacities = Array.from(triangles).map((t) => parseFloat(t.getAttribute("fill-opacity") ?? "1"));
    expect(opacities).toContain(0.3);
    // The non-excluded one stays bright (>= 0.9)
    expect(Math.max(...opacities)).toBeGreaterThanOrEqual(0.9);
  });
});

describe("<TraceViewer> — overlay ticks", () => {
  it("renders one overlay <line> per predicted_q in the active group", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const indices: IndexEntry[] = [{
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7, 0.9],
      peaks: [],
    }];
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        activeGroupIndices={indices}
        hoveredIndex={undefined}
        {...defaultProps}
      />,
    );
    const ticks = container.querySelectorAll('[data-role="tick-root"] line');
    expect(ticks.length).toBe(2);
    // No detected peaks → ticks are unmatched → dimmed to 0.15 so the user
    // can see which predicted positions were not confirmed by the data.
    const first = ticks[0] as SVGLineElement;
    expect(first.getAttribute("stroke-opacity")).toBe("0.15");
  });

  it("matched ticks (detected peak at same q) render at full base opacity", () => {
    const trace = { q: [0.1, 0.5, 0.9], I: [10, 20, 10], sigma: [1, 1, 1] };
    const indices: IndexEntry[] = [{
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.5],
      peaks: [],
    }];
    // Peak at exactly q=0.5 — same pixel as the predicted_q → matched.
    const matchedPeak = {
      id: 99, exposure_id: 1, q: 0.5, intensity: 20, prominence: null,
      sharpness: null, source: "auto" as const, excluded: false,
    };
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[matchedPeak]}
        activeGroupIndices={indices}
        hoveredIndex={undefined}
        {...defaultProps}
      />,
    );
    const ticks = container.querySelectorAll('[data-role="tick-root"] line');
    expect(ticks.length).toBe(1);
    expect((ticks[0] as SVGLineElement).getAttribute("stroke-opacity")).toBe("0.35");
  });

  it("adds the hovered index's ticks on top of active ones", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
    const active: IndexEntry[] = [{
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7, 0.9],
      peaks: [],
    }];
    const hovered: IndexEntry = {
      id: 11, exposure_id: 1, phase: "Im3m", basis: 0.3, score: 0.5,
      r_squared: 0.7, lattice_d: 9.0, status: "candidate",
      predicted_q: [0.42, 0.6, 0.8],
      peaks: [],
    };
    const { container } = render(
      <TraceViewer
        trace={trace}
        peaks={[]}
        activeGroupIndices={active}
        hoveredIndex={hovered}
        {...defaultProps}
      />,
    );
    const ticks = container.querySelectorAll('[data-role="tick-root"] line');
    // 2 active + 3 hovered = 5 ticks total
    expect(ticks.length).toBe(5);
    // Active (faded, unmatched): neutral gray at 0.12.
    const dimmed = Array.from(ticks).slice(0, 2) as SVGLineElement[];
    const dimmedOpacities = dimmed.map((n) => n.getAttribute("stroke-opacity"));
    const dimmedStrokes   = dimmed.map((n) => n.getAttribute("stroke"));
    expect(dimmedOpacities.every((o) => o === "0.12")).toBe(true);
    expect(dimmedStrokes.every((s) => s === "var(--color-fg-dim)")).toBe(true);
    // Hovered (strong, unmatched): phase color at 0.45.
    const hov = Array.from(ticks).slice(2) as SVGLineElement[];
    expect(hov.every((n) => n.getAttribute("stroke-opacity") === "0.45")).toBe(true);
  });
});
