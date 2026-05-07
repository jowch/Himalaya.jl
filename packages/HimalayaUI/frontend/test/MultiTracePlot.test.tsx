/**
 * MultiTracePlot tests (Plan §Phase 6, Task 6.2).
 *
 * Two layers of testing:
 *  1. Pure y-band layout math (`computeYBands`) — runs in JSDOM trivially.
 *  2. Component rendering — Observable Plot is mocked at the module level
 *     (mirrors `MillerPlot.test.tsx`). The mock captures the marks array
 *     plus event handler registration.
 *
 * Tests pin behavior — they don't pin styling. Brush/double-click semantics
 * are exercised by simulating the registered event listeners on the mock
 * plot element (the same pattern Observable Plot uses internally).
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render } from "@testing-library/react";
import { MultiTracePlot, computeYBands, COMPARE_PLOT_ASPECT } from "../src/components/MultiTracePlot";
import type { ComparisonMember } from "../src/api";

let lastPlotElement: HTMLElement | null = null;

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "multi-trace-svg");
    // Stub the runtime `.scale("x")` interface used by invertQ.
    (el as unknown as { scale: (n: string) => { invert?: (px: number) => number; apply?: (q: number) => number } | undefined }).scale = (n) =>
      n === "x"
        ? { invert: (px: number) => px / 100, apply: (q: number) => q * 100 }
        : undefined;
    lastPlotElement = el;
    return el;
  }),
  line: vi.fn((data: unknown, opts: unknown) => ({ _kind: "line", data, opts })),
  dot:  vi.fn((data: unknown, opts: unknown) => ({ _kind: "dot",  data, opts })),
  text: vi.fn((data: unknown, opts: unknown) => ({ _kind: "text", data, opts })),
}));

beforeEach(() => {
  lastPlotElement = null;
});

// ── computeYBands ────────────────────────────────────────────────────────────

describe("computeYBands", () => {
  it("returns one band per member, summing to panelHeight", () => {
    const bands = computeYBands([1, 1, 1], 300);
    expect(bands).toHaveLength(3);
    expect(bands[0]).toEqual([0, 100]);
    expect(bands[1]).toEqual([100, 200]);
    expect(bands[2]).toEqual([200, 300]);
  });

  it("respects unequal band_height ratios", () => {
    const bands = computeYBands([1, 2, 1], 400);
    // ratio sum = 4; each ratio*100 = pixel slice.
    expect(bands[0]).toEqual([0, 100]);
    expect(bands[1]).toEqual([100, 300]);
    expect(bands[2]).toEqual([300, 400]);
  });

  it("reorders bands top-down according to caller-supplied order (display_order)", () => {
    // Same ratios but if display_order changes, the caller passes the
    // reordered ratios; the function preserves input order.
    const a = computeYBands([1, 2], 300);
    const b = computeYBands([2, 1], 300);
    expect(a[0]![1]).toBe(100);  // first band 100 px tall when ratio = 1
    expect(b[0]![1]).toBe(200);  // first band 200 px tall when ratio = 2
  });

  it("returns an empty array for no members", () => {
    expect(computeYBands([], 300)).toEqual([]);
  });

  it("handles a zero panel height gracefully (degenerate render frame)", () => {
    const bands = computeYBands([1, 1], 0);
    expect(bands[0]).toEqual([0, 0]);
    expect(bands[1]).toEqual([0, 0]);
  });
});

// ── COMPARE_PLOT_ASPECT ──────────────────────────────────────────────────────

describe("COMPARE_PLOT_ASPECT", () => {
  it("is hardcoded at 0.3 per spec §Plot rendering", () => {
    expect(COMPARE_PLOT_ASPECT).toBe(0.3);
  });
});

// ── MultiTracePlot ───────────────────────────────────────────────────────────

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1, comparison_id: 100, exposure_id: 42, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: null, created_at: null,
    ...over,
  };
}

describe("<MultiTracePlot>", () => {
  it("renders a container with testid 'multi-trace-plot'", () => {
    const { container } = render(
      <MultiTracePlot
        members={[]}
        traces={new Map()}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );
    expect(container.querySelector('[data-testid="multi-trace-plot"]')).not.toBeNull();
  });

  it("renders an empty plot for zero members (does not crash)", () => {
    expect(() =>
      render(
        <MultiTracePlot
          members={[]}
          traces={new Map()}
          xDomain={null}
          onXDomain={() => {}}
        />,
      ),
    ).not.toThrow();
  });

  it("renders a Plot with one mark group per member when members + traces are present", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.plot as unknown as { mockClear: () => void }).mockClear();
    (Plot.line as unknown as { mockClear: () => void }).mockClear();

    const m1 = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const m2 = makeMember({ id: 2, exposure_id: 20, display_order: 1 });
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [0, 0, 0] };
    const traces = new Map([[10, trace], [20, trace]]);

    render(
      <MultiTracePlot
        members={[m1, m2]}
        traces={traces}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    expect(Plot.plot).toHaveBeenCalledTimes(1);
    // Two members, each contributes a line mark.
    expect(Plot.line).toHaveBeenCalledTimes(2);
  });

  it("calls onXDomain(null) on double-click (reset zoom)", async () => {
    const onXDomain = vi.fn();
    const m1 = makeMember({ id: 1, exposure_id: 10 });
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [0, 0] };

    render(
      <MultiTracePlot
        members={[m1]}
        traces={new Map([[10, trace]])}
        xDomain={[0.1, 0.2]}
        onXDomain={onXDomain}
      />,
    );

    // Trigger a dblclick on the captured plot element.
    expect(lastPlotElement).not.toBeNull();
    const ev = new MouseEvent("dblclick", { bubbles: true });
    lastPlotElement!.dispatchEvent(ev);

    expect(onXDomain).toHaveBeenCalledWith(null);
  });

  it("re-renders the plot when members reorder (regression: bands shift)", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.plot as unknown as { mockClear: () => void }).mockClear();

    const m1 = makeMember({ id: 1, exposure_id: 10, display_order: 0, band_height: 1 });
    const m2 = makeMember({ id: 2, exposure_id: 20, display_order: 1, band_height: 1 });
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [0, 0] };
    const traces = new Map([[10, trace], [20, trace]]);

    const { rerender } = render(
      <MultiTracePlot
        members={[m1, m2]}
        traces={traces}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    // Reorder: swap display_order so m2 is on top.
    const m1Reordered = { ...m1, display_order: 1 };
    const m2Reordered = { ...m2, display_order: 0 };
    rerender(
      <MultiTracePlot
        members={[m2Reordered, m1Reordered]}
        traces={traces}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    // Plot.plot should have been re-invoked because props changed.
    expect((Plot.plot as unknown as { mock: { calls: unknown[][] } }).mock.calls.length).toBeGreaterThanOrEqual(2);
  });
});
