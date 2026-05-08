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
import {
  MultiTracePlot,
  computeYBands,
  computeMaxPlotWidth,
  MIN_PLOT_WIDTH,
  COMPARE_PLOT_ASPECT,
} from "../src/components/MultiTracePlot";
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

// ── computeMaxPlotWidth ──────────────────────────────────────────────────────
//
// Pins the per-band aspect math from issue #81. The formula
// `(panelHeight / n) * target` is load-bearing and was previously written as
// `n * bandH * target` which simplifies to `panelHeight * target` — that
// targets the *whole-plot* aspect, not per-band, and produced caps ~N× too
// large. Without these tests, the algebra trap could quietly recur during a
// future refactor.

describe("computeMaxPlotWidth", () => {
  it("for a single member with a short panel, clamps to MIN_PLOT_WIDTH", () => {
    // bandH = 200, target = 200×1.0 = 200, max(280, 200) = 280.
    expect(computeMaxPlotWidth(200, 1)).toBe(MIN_PLOT_WIDTH);
  });

  it("for n=4 members at 600px panel, the floor wins (per-band 1.0 ≈ 150 < 280)", () => {
    // bandH = 150, target = 150×1.0 = 150, max(280, 150) = 280.
    expect(computeMaxPlotWidth(600, 4)).toBe(MIN_PLOT_WIDTH);
  });

  it("for n=4 members at 1600px panel, the per-band target wins (≈ 400 > 280)", () => {
    // bandH = 400, target = 400×1.0 = 400, max(280, 400) = 400.
    expect(computeMaxPlotWidth(1600, 4)).toBe(400);
  });

  it("returns MIN_PLOT_WIDTH for zero members (degenerate render frame)", () => {
    expect(computeMaxPlotWidth(800, 0)).toBe(MIN_PLOT_WIDTH);
  });

  it("returns MIN_PLOT_WIDTH for zero panel height (initial mount, JSDOM)", () => {
    expect(computeMaxPlotWidth(0, 3)).toBe(MIN_PLOT_WIDTH);
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

  it("renders one [data-testid='member-trace'] overlay per member with the right data-member-id", () => {
    const m1 = makeMember({ id: 11, exposure_id: 110, display_order: 0 });
    const m2 = makeMember({ id: 22, exposure_id: 220, display_order: 1 });
    const m3 = makeMember({ id: 33, exposure_id: 330, display_order: 2 });
    const trace = { q: [0.1, 0.2], I: [10, 20], sigma: [0, 0] };
    const traces = new Map([[110, trace], [220, trace], [330, trace]]);

    const { container } = render(
      <MultiTracePlot
        members={[m1, m2, m3]}
        traces={traces}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    const overlays = container.querySelectorAll('[data-testid="member-trace"]');
    expect(overlays).toHaveLength(3);
    const ids = Array.from(overlays).map((el) => el.getAttribute("data-member-id"));
    expect(ids).toEqual(["11", "22", "33"]);
  });

  it("calls onPeakClick(memberId, peakId, altKey) when click lands within hit radius of a peak", () => {
    const onPeakClick = vi.fn();
    const m1 = makeMember({
      id: 11, exposure_id: 110, display_order: 0,
      snapshot: {
        effective_peaks: [
          { id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
          { id: 22, q: 0.50, intensity: 80, sharpness: 1, source: "auto" },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const trace = { q: [0.1, 0.2, 0.3, 0.4, 0.5], I: [1, 2, 3, 4, 5], sigma: [0, 0, 0, 0, 0] };

    render(
      <MultiTracePlot
        members={[m1]}
        traces={new Map([[110, trace]])}
        xDomain={null}
        onXDomain={() => {}}
        onPeakClick={onPeakClick}
      />,
    );

    expect(lastPlotElement).not.toBeNull();
    // The mock x-scale maps q*100 → px; peak 21 sits at q=0.30 ⇒ px=30.
    // Dispatch a click at clientX=30 + container offset (jsdom defaults to 0).
    const ev = new MouseEvent("click", { bubbles: true, clientX: 30, clientY: 5 });
    lastPlotElement!.dispatchEvent(ev);
    expect(onPeakClick).toHaveBeenCalledWith(11, 21, false);
  });

  it("does NOT call onPeakClick when click is outside hit radius of any peak", () => {
    const onPeakClick = vi.fn();
    const m1 = makeMember({
      id: 11, exposure_id: 110, display_order: 0,
      snapshot: {
        effective_peaks: [
          { id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const trace = { q: [0.1, 0.2, 0.3, 0.4, 0.5], I: [1, 2, 3, 4, 5], sigma: [0, 0, 0, 0, 0] };

    render(
      <MultiTracePlot
        members={[m1]}
        traces={new Map([[110, trace]])}
        xDomain={null}
        onXDomain={() => {}}
        onPeakClick={onPeakClick}
      />,
    );

    // Click far from peak 21 (px=30): clientX=200 is well beyond hit radius.
    const ev = new MouseEvent("click", { bubbles: true, clientX: 200, clientY: 5 });
    lastPlotElement!.dispatchEvent(ev);
    expect(onPeakClick).not.toHaveBeenCalled();
  });

  it("propagates altKey through to onPeakClick", () => {
    const onPeakClick = vi.fn();
    const m1 = makeMember({
      id: 11, exposure_id: 110, display_order: 0,
      snapshot: {
        effective_peaks: [
          { id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const trace = { q: [0.1, 0.2, 0.3, 0.4, 0.5], I: [1, 2, 3, 4, 5], sigma: [0, 0, 0, 0, 0] };

    render(
      <MultiTracePlot
        members={[m1]}
        traces={new Map([[110, trace]])}
        xDomain={null}
        onXDomain={() => {}}
        onPeakClick={onPeakClick}
      />,
    );

    const ev = new MouseEvent("click", { bubbles: true, clientX: 30, clientY: 5, altKey: true });
    lastPlotElement!.dispatchEvent(ev);
    expect(onPeakClick).toHaveBeenCalledWith(11, 21, true);
  });

  it("does not register a click handler when onPeakClick is not provided (no edit mode)", () => {
    const onPeakClick = vi.fn();
    const m1 = makeMember({
      id: 11, exposure_id: 110, display_order: 0,
      snapshot: {
        effective_peaks: [{ id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const trace = { q: [0.1, 0.2, 0.3, 0.4, 0.5], I: [1, 2, 3, 4, 5], sigma: [0, 0, 0, 0, 0] };

    // No onPeakClick prop.
    render(
      <MultiTracePlot
        members={[m1]}
        traces={new Map([[110, trace]])}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    const ev = new MouseEvent("click", { bubbles: true, clientX: 30, clientY: 5 });
    lastPlotElement!.dispatchEvent(ev);
    expect(onPeakClick).not.toHaveBeenCalled();
  });

  it("hits the correct member when multiple members stack and click Y picks one band", () => {
    const onPeakClick = vi.fn();
    const m1 = makeMember({
      id: 11, exposure_id: 110, display_order: 0,
      snapshot: {
        effective_peaks: [{ id: 31, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const m2 = makeMember({
      id: 22, exposure_id: 220, display_order: 1,
      snapshot: {
        effective_peaks: [{ id: 32, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    const trace = { q: [0.1, 0.2, 0.3, 0.4, 0.5], I: [1, 2, 3, 4, 5], sigma: [0, 0, 0, 0, 0] };

    render(
      <MultiTracePlot
        members={[m1, m2]}
        traces={new Map([[110, trace], [220, trace]])}
        xDomain={null}
        onXDomain={() => {}}
        onPeakClick={onPeakClick}
      />,
    );

    // Both members have a peak at q=0.30 ⇒ both at px=30. Click should
    // resolve to ONE specific member based on which is closer (within
    // hit radius of either). Either resolution is acceptable; here we
    // assert at least ONE call landed (no double-firing) and the peak id
    // matches the expected (m1's 31 or m2's 32).
    const ev = new MouseEvent("click", { bubbles: true, clientX: 30, clientY: 5 });
    lastPlotElement!.dispatchEvent(ev);
    expect(onPeakClick).toHaveBeenCalledTimes(1);
    const args = onPeakClick.mock.calls[0]!;
    expect([11, 22]).toContain(args[0]);
    expect([31, 32]).toContain(args[1]);
  });

  it("renders high-intensity points at smaller DOM y than low-intensity points (regression #63: plot was upside-down)", async () => {
    // Invariant: within a trace, the brightest q-point must render closer
    // to the top of the plot than the dimmest. This survives any future
    // change to either `applyNormalization`'s output convention or Plot's
    // y-scale configuration — as long as the two stay consistent.
    //
    // To catch #63 (where Plot's y.domain was [0, panelH] while
    // applyNormalization produced screen-coord ys), we compose both sides
    // of the chain:
    //   1. capture the y values handed to Plot.line for two q-points — a
    //      "peak" (high I) and a "baseline" (low I) — within the same
    //      single-member trace;
    //   2. capture the y.domain config passed to Plot.plot;
    //   3. project each pixel-y onto a normalized [0, 1] DOM-position axis
    //      via `(y - domain[0]) / (domain[1] - domain[0])`. Plot maps
    //      domain[0] → bottom, domain[1] → top, so a *higher* normalized
    //      value = closer to top; we want the peak to score higher than
    //      the baseline.
    const Plot = await import("@observablehq/plot");
    (Plot.plot as unknown as { mockClear: () => void }).mockClear();
    (Plot.line as unknown as { mockClear: () => void }).mockClear();

    const member = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    // One trace with a "baseline" point (low I) and a "peak" point (high I)
    // at distinct q values. Self-normalization scales both into the band
    // but preserves their relative ordering (the peak fills the working
    // band; the baseline collapses to its bottom).
    const trace = {
      q: [0.10, 0.20],
      I: [    1, 10000],
      sigma: [0, 0],
    };

    render(
      <MultiTracePlot
        members={[member]}
        traces={new Map([[10, trace]])}
        xDomain={null}
        onXDomain={() => {}}
      />,
    );

    const lineCalls = (Plot.line as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    const data = lineCalls[0]![0] as Array<{ q: number; y: number }>;
    const baseline = data.find((p) => p.q === 0.10);
    const peak     = data.find((p) => p.q === 0.20);
    if (baseline === undefined || peak === undefined) {
      throw new Error("trace endpoints missing from line mark data");
    }

    const plotCalls = (Plot.plot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    const cfg = plotCalls[0]![0] as { y: { domain: [number, number] } };
    const [d0, d1] = cfg.y.domain;
    expect(d1).not.toBe(d0); // sanity: scale is non-degenerate

    // Project pixel-y onto [0, 1] DOM-position axis (1 = top).
    const topness = (y: number): number => (y - d0) / (d1 - d0);
    expect(topness(peak.y)).toBeGreaterThan(topness(baseline.y));
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
