/**
 * AnnotationToggles tests (Plan §Phase 9, Task 9.3).
 *
 * - Both toggles render with correct default state (true / true).
 * - Toggling fires the named Zustand action.
 * - When the global flag is off, `buildMemberMarks` doesn't emit triangles
 *   (peaks toggle) or text labels (labels toggle).
 * - The toggle is mounted by ComparePage (review) but not by ComparePageEdit
 *   (edit-mode peak click cycle handles per-peak control instead).
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as Plot from "@observablehq/plot";
import { render, screen, fireEvent } from "@testing-library/react";
import { AnnotationToggles } from "../src/components/AnnotationToggles";
import { useAppState } from "../src/state";
import { buildMemberMarks } from "../src/components/MemberTraceLayer";
import type { ComparisonMember } from "../src/api";

vi.mock("@observablehq/plot", () => ({
  // Stub Plot.plot for the ComparePage mount tests (MultiTracePlot calls it).
  plot: vi.fn(() => {
    const el = document.createElement("div");
    (el as unknown as { scale: (n: string) => unknown }).scale = (n) =>
      n === "x"
        ? { invert: (px: number) => px / 100, apply: (q: number) => q * 100 }
        : undefined;
    return el;
  }),
  line: vi.fn((data: unknown, opts: unknown) => ({ _kind: "line", data, opts })),
  dot:  vi.fn((data: unknown, opts: unknown) => ({ _kind: "dot",  data, opts })),
  text: vi.fn((data: unknown, opts: unknown) => ({ _kind: "text", data, opts })),
  link: vi.fn((data: unknown, opts: unknown) => ({ _kind: "link", data, opts })),
}));

beforeEach(() => {
  useAppState.setState({ showPeakTicks: true, showPeakLabels: true });
  (Plot.line as unknown as { mockClear: () => void }).mockClear();
  (Plot.dot  as unknown as { mockClear: () => void }).mockClear();
  (Plot.text as unknown as { mockClear: () => void }).mockClear();
});

describe("AnnotationToggles — render + Zustand wiring", () => {
  it("both toggle buttons render with default state active=true", () => {
    render(<AnnotationToggles />);
    const peaks = screen.getByTestId("annotation-toggle-peaks");
    const labels = screen.getByTestId("annotation-toggle-labels");
    expect(peaks).toHaveAttribute("data-active", "true");
    expect(peaks).toHaveAttribute("aria-pressed", "true");
    expect(labels).toHaveAttribute("data-active", "true");
    expect(labels).toHaveAttribute("aria-pressed", "true");
  });

  it("toggling peaks fires setShowPeakTicks", () => {
    render(<AnnotationToggles />);
    const peaks = screen.getByTestId("annotation-toggle-peaks");
    fireEvent.click(peaks);
    expect(useAppState.getState().showPeakTicks).toBe(false);
  });

  it("toggling labels fires setShowPeakLabels", () => {
    render(<AnnotationToggles />);
    const labels = screen.getByTestId("annotation-toggle-labels");
    fireEvent.click(labels);
    expect(useAppState.getState().showPeakLabels).toBe(false);
  });

  it("toggle state reflects the store after external set", () => {
    useAppState.setState({ showPeakTicks: false });
    render(<AnnotationToggles />);
    const peaks = screen.getByTestId("annotation-toggle-peaks");
    expect(peaks).toHaveAttribute("data-active", "false");
    expect(peaks).toHaveAttribute("aria-pressed", "false");
  });
});

// Mark-level integration: the toggles control mark emission via
// `buildMemberMarks`, not the rendered DOM (Observable Plot is mocked).
function makeMember(): ComparisonMember {
  return {
    id: 1, comparison_id: 100, exposure_id: 42,
    display_order: 0, band_height: 1, y_offset: 0,
    normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null,
    peak_display: { hidden: [], labeled: [11, 12] },
    snapshot: {
      effective_peaks: [
        { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        { id: 12, q: 0.50, intensity: 80, sharpness: 1, source: "auto" },
      ],
      confirmed_index: {
        id: 7, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.5,
        peak_ids: [11, 12],
      },
      analysis_inputs_hash: "abc",
    },
    is_stale: false, created_by: null, created_at: null,
  };
}
const trace = {
  q: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
  I: [10, 20, 30, 50, 80, 60],
  sigma: [0, 0, 0, 0, 0, 0],
};

describe("AnnotationToggles — buildMemberMarks honors flags", () => {
  it("showPeakTicks=false suppresses dot marks", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12] },
      showPeakTicks: false,
      showPeakLabels: true,
    });
    expect(Plot.dot).not.toHaveBeenCalled();
    // Line is unaffected.
    expect(Plot.line).toHaveBeenCalledTimes(1);
  });

  it("showPeakLabels=false suppresses text marks but keeps dots", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12] },
      showPeakTicks: true,
      showPeakLabels: false,
    });
    // Dots still rendered (peaks toggle on)
    expect((Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls.length)
      .toBeGreaterThanOrEqual(1);
    // No text labels (labels toggle off)
    expect(Plot.text).not.toHaveBeenCalled();
  });

  it("both off → only the line mark", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12] },
      showPeakTicks: false,
      showPeakLabels: false,
    });
    expect(Plot.line).toHaveBeenCalledTimes(1);
    expect(Plot.dot).not.toHaveBeenCalled();
    expect(Plot.text).not.toHaveBeenCalled();
  });

  it("defaults (both true / omitted) emit dots + labels normally", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12] },
    });
    expect((Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls.length)
      .toBeGreaterThanOrEqual(1);
    expect((Plot.text as unknown as { mock: { calls: unknown[][] } }).mock.calls.length)
      .toBeGreaterThanOrEqual(1);
  });
});

describe("AnnotationToggles — page mount surface", () => {
  it("ComparePage (review mode) mounts the toggles", async () => {
    const { MemoryRouter, Routes, Route } = await import("react-router-dom");
    const { QueryClient, QueryClientProvider } = await import("@tanstack/react-query");
    const { ComparePage } = await import("../src/pages/ComparePage");
    const qc = new QueryClient({
      defaultOptions: {
        queries: { retry: false, gcTime: Infinity, staleTime: 0 },
        mutations: { retry: false },
      },
    });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("annotation-toggles")).toBeInTheDocument();
  });

  it("ComparePageEdit (edit mode) does NOT mount the toggles", async () => {
    const { MemoryRouter, Routes, Route } = await import("react-router-dom");
    const { QueryClient, QueryClientProvider } = await import("@tanstack/react-query");
    const { ComparePageEdit } = await import("../src/pages/ComparePageEdit");
    const qc = new QueryClient({
      defaultOptions: {
        queries: { retry: false, gcTime: Infinity, staleTime: 0 },
        mutations: { retry: false },
      },
    });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/new"]}>
          <Routes>
            <Route path="/experiments/:eid/compare/new" element={<ComparePageEdit />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.queryByTestId("annotation-toggles")).toBeNull();
  });
});
