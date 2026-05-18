/**
 * Peak hover tooltip tests (Plan §Phase 8, Task 8.3).
 *
 * Triangles get a tooltip showing q-value on hover. When a developer-only
 * `?showPeakIds` URL flag is set, the tooltip also includes the peak id
 * (kept out of production UI to avoid leaking implementation detail).
 *
 * Test surface: simulate mousemove + mouseleave on the captured plot
 * element; assert against a `data-testid="peak-tooltip"` overlay rendered
 * by `MultiTracePlot`. Hit-testing uses the same pixel-radius logic as
 * the click handler (Phase 8.1).
 */
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { act, render, screen, waitFor } from "@testing-library/react";
import { MultiTracePlot } from "../src/components/MultiTracePlot";
import type { SeriesMember } from "../src/api";

let lastPlotElement: HTMLElement | null = null;

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "multi-trace-svg");
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
  link: vi.fn((data: unknown, opts: unknown) => ({ _kind: "link", data, opts })),
}));

beforeEach(() => {
  lastPlotElement = null;
});

afterEach(() => {
  // Reset the URL-search state so flag tests don't leak.
  if (typeof window !== "undefined") {
    window.history.replaceState({}, "", "/");
  }
});

function makeMember(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1, series_id: 100, exposure_id: 42, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: null, created_at: null,
    ...over,
  };
}

describe("peak hover tooltip (Phase 8.3)", () => {
  it("shows tooltip with q-value when mouse hovers near a peak triangle", async () => {
    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
        onPeakClick={() => {}}
      />,
    );

    // Hover at q=0.30 → px=30 (mock scale apply = q*100).
    expect(lastPlotElement).not.toBeNull();
    const ev = new MouseEvent("mousemove", { bubbles: true, clientX: 30, clientY: 5 });
    lastPlotElement!.dispatchEvent(ev);

    const tip = await screen.findByTestId("peak-tooltip");
    expect(tip).not.toBeNull();
    // q rendered to 3 sig digits.
    expect(tip.textContent).toMatch(/0\.300?/);
  });

  it("hides tooltip on mouseleave", async () => {
    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
        onPeakClick={() => {}}
      />,
    );

    act(() => {
      lastPlotElement!.dispatchEvent(new MouseEvent("mousemove", { bubbles: true, clientX: 30, clientY: 5 }));
    });
    await screen.findByTestId("peak-tooltip");

    act(() => {
      lastPlotElement!.dispatchEvent(new MouseEvent("mouseleave", { bubbles: true }));
    });
    await waitFor(() => {
      expect(screen.queryByTestId("peak-tooltip")).toBeNull();
    });
  });

  it("does NOT show tooltip when mouse is far from any peak", () => {
    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
        onPeakClick={() => {}}
      />,
    );

    lastPlotElement!.dispatchEvent(new MouseEvent("mousemove", { bubbles: true, clientX: 200, clientY: 5 }));
    expect(screen.queryByTestId("peak-tooltip")).toBeNull();
  });

  it("includes peak id when ?showPeakIds URL flag is set", async () => {
    window.history.replaceState({}, "", "/?showPeakIds");

    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 99, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
        onPeakClick={() => {}}
      />,
    );

    act(() => {
      lastPlotElement!.dispatchEvent(new MouseEvent("mousemove", { bubbles: true, clientX: 30, clientY: 5 }));
    });
    const tip = await screen.findByTestId("peak-tooltip");
    expect(tip.textContent).toContain("99");
  });

  it("omits peak id when the URL flag is NOT set (production UI)", async () => {
    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 99, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
        onPeakClick={() => {}}
      />,
    );

    act(() => {
      lastPlotElement!.dispatchEvent(new MouseEvent("mousemove", { bubbles: true, clientX: 30, clientY: 5 }));
    });
    const tip = await screen.findByTestId("peak-tooltip");
    expect(tip.textContent).not.toContain("99");
  });

  it("renders tooltip even in review mode (no onPeakClick) — tooltip is unconditional", async () => {
    const m1 = makeMember({
      id: 11, exposure_id: 110,
      snapshot: {
        effective_peaks: [{ id: 21, q: 0.30, intensity: 50, sharpness: 1, source: "auto" }],
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
      />,
    );

    act(() => {
      lastPlotElement!.dispatchEvent(new MouseEvent("mousemove", { bubbles: true, clientX: 30, clientY: 5 }));
    });
    const tip = await screen.findByTestId("peak-tooltip");
    expect(tip.textContent).toMatch(/0\.300?/);
  });
});
