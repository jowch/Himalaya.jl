import { describe, it, expect, vi } from "vitest";
import { render, fireEvent, cleanup } from "@testing-library/react";
import { afterEach, beforeAll, afterAll } from "vitest";
import { PlotSurface } from "../src/components/PlotSurface";
import type { Peak } from "../src/api";

// Mock Observable Plot: Plot.plot returns an element carrying a stub `.scale`
// so the surface's scale plumbing (apply/invert) round-trips deterministically.
// apply: q → q*1000 px; invert: px → px/1000 q. Matches TraceViewer's mock.
vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "plot-svg");
    // getBoundingClientRect drives the gesture coordinate math.
    el.getBoundingClientRect = () =>
      ({ left: 0, top: 0, width: 400, height: 300, right: 400, bottom: 300, x: 0, y: 0, toJSON() {} }) as DOMRect;
    (el as unknown as { scale: (n: string) => unknown }).scale = (name: string) =>
      name === "x"
        ? { apply: (q: number) => q * 1000, invert: (px: number) => px / 1000 }
        : { apply: (v: number) => v, invert: (px: number) => px };
    return el;
  }),
  line: vi.fn(() => ({ _kind: "line" })),
  areaY: vi.fn(() => ({ _kind: "areaY" })),
  dot: vi.fn(() => ({ _kind: "dot" })),
}));

const RECT = { left: 0, top: 0, width: 400, height: 300, right: 400, bottom: 300, x: 0, y: 0, toJSON() {} } as DOMRect;
let origRect: typeof Element.prototype.getBoundingClientRect;
beforeAll(() => {
  origRect = Element.prototype.getBoundingClientRect;
  Element.prototype.getBoundingClientRect = () => RECT;
});
afterAll(() => {
  Element.prototype.getBoundingClientRect = origRect;
});
afterEach(() => cleanup());

const MARGINS = { left: 50, right: 14, top: 36, bottom: 40 };

function peak(id: number, q: number): Peak {
  return {
    id, exposure_id: 1, q, intensity: null, prominence: null,
    sharpness: null, source: "auto", excluded: false,
  };
}

describe("<PlotSurface>", () => {
  it("renders the data marks via Plot.plot", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.plot as ReturnType<typeof vi.fn>).mockClear();
    render(
      <PlotSurface
        marks={[{ _kind: "line" } as never]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        onXDomain={() => {}}
      />,
    );
    expect(Plot.plot).toHaveBeenCalled();
    const opts = (Plot.plot as ReturnType<typeof vi.fn>).mock.calls.at(-1)![0] as {
      marginLeft: number; x: { type: string };
    };
    expect(opts.marginLeft).toBe(50);
    expect(opts.x.type).toBe("log");
  });

  it("a wheel event narrows the x-domain via onXDomain", () => {
    const onXDomain = vi.fn();
    const { container } = render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        onXDomain={onXDomain}
      />,
    );
    const plotEl = container.querySelector('[data-testid="plot-svg"]')!;
    fireEvent.wheel(plotEl, { clientX: 200, clientY: 100, deltaY: -100 });
    expect(onXDomain).toHaveBeenCalled();
    const dom = onXDomain.mock.calls.at(-1)![0] as [number, number];
    expect(dom[0]).toBeGreaterThanOrEqual(0.01);
    expect(dom[1]).toBeLessThanOrEqual(1);
    expect(dom[1] - dom[0]).toBeLessThan(1 - 0.01); // narrowed
  });

  it("double-click resets the domain (onXDomain(null) by default)", () => {
    const onXDomain = vi.fn();
    const { container } = render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.1, 0.5]}
        margins={MARGINS}
        onXDomain={onXDomain}
      />,
    );
    const plotEl = container.querySelector('[data-testid="plot-svg"]')!;
    fireEvent.dblClick(plotEl);
    expect(onXDomain).toHaveBeenLastCalledWith(null);
  });

  it("a click within tolerance of a peak calls onClickPeak", () => {
    const onClickPeak = vi.fn();
    const onAddPeak = vi.fn();
    const { container } = render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        peaks={[peak(7, 0.1)]} // px = 100
        onXDomain={() => {}}
        onClickPeak={onClickPeak}
        onAddPeak={onAddPeak}
      />,
    );
    const plotEl = container.querySelector('[data-testid="plot-svg"]')!;
    // click at px=102 (within PEAK_HIT_PX of the peak at 100), y inside interior
    fireEvent.click(plotEl, { clientX: 102, clientY: 150 });
    expect(onClickPeak).toHaveBeenCalledWith(7, false);
    expect(onAddPeak).not.toHaveBeenCalled();
  });

  it("an empty-area click calls onAddPeak with the inverted q", () => {
    const onClickPeak = vi.fn();
    const onAddPeak = vi.fn();
    const { container } = render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        peaks={[peak(7, 0.1)]}
        onXDomain={() => {}}
        onClickPeak={onClickPeak}
        onAddPeak={onAddPeak}
      />,
    );
    const plotEl = container.querySelector('[data-testid="plot-svg"]')!;
    // click at px=300 (far from the peak), interior y → add a peak at q=0.3
    fireEvent.click(plotEl, { clientX: 300, clientY: 150 });
    expect(onAddPeak).toHaveBeenCalled();
    expect(onAddPeak.mock.calls.at(-1)![0]).toBeCloseTo(0.3, 5);
    expect(onClickPeak).not.toHaveBeenCalled();
  });

  it("never routes a click to an optimistic (id<0) peak", () => {
    const onClickPeak = vi.fn();
    const onAddPeak = vi.fn();
    const { container } = render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        peaks={[peak(-3, 0.1)]} // optimistic placeholder at px=100
        onXDomain={() => {}}
        onClickPeak={onClickPeak}
        onAddPeak={onAddPeak}
      />,
    );
    const plotEl = container.querySelector('[data-testid="plot-svg"]')!;
    fireEvent.click(plotEl, { clientX: 100, clientY: 150 });
    expect(onClickPeak).not.toHaveBeenCalled();
    // it falls through to an add at the clicked q
    expect(onAddPeak).toHaveBeenCalled();
  });

  it("overlay(ctx) receives applyQ/invertQ + a hitTest that skips optimistic peaks", () => {
    let captured: import("../src/components/PlotSurface").PlotOverlayContext | null = null;
    render(
      <PlotSurface
        marks={[]}
        xType="log"
        xDomain={[0.01, 1]}
        margins={MARGINS}
        onXDomain={() => {}}
        overlay={(ctx) => {
          captured = ctx;
          return null;
        }}
      />,
    );
    expect(captured).not.toBeNull();
    const ctx = captured!;
    expect(ctx.applyQ(0.1)).toBeCloseTo(100, 5);
    expect(ctx.invertQ(100)).toBeCloseTo(0.1, 5);
    // hitTest finds the real peak…
    expect(ctx.hitTest([peak(7, 0.1)], 100, 10)?.id).toBe(7);
    // …and skips the optimistic one.
    expect(ctx.hitTest([peak(-3, 0.1)], 100, 10)).toBeNull();
  });
});
