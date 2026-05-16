/**
 * Compare page shells (Plan §Phase 4, Task 4.1; merged by Compare UX C-15).
 *
 * Verifies that the unified `Compare` component renders the review body and
 * the edit body under the expected routes and reads URL params correctly.
 * `Compare` picks the edit body when the URL ends `/new` or a draft is active,
 * and the review body otherwise — so the edit-mode assertions seed a draft.
 *
 * Behaviour exercised by later tasks (sidebar, draft state, save flow) lives
 * in their own test files.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor } from "@testing-library/react";
import { Compare } from "../src/pages/Compare";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";

function renderAt(path: string) {
  const qc = new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/experiments/:eid/compare" element={<Compare />} />
          <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          <Route path="/experiments/:eid/compare/new" element={<Compare />} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<Compare />} />
          <Route path="/compare/all" element={<Compare />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

/** A minimal draft tied to `id` so `Compare` mounts the edit body. */
function seedDraftFor(id: number | undefined): void {
  useAppState.setState({
    activeDraft: {
      id,
      baseHash: id !== undefined ? "h" : undefined,
      title: "T",
      description: "",
      members: [],
      forkedFromId: undefined,
      forkedAtHash: undefined,
      viewGroupingMode: undefined,
      viewShowPeakTicks: undefined,
      viewShowPeakLabels: undefined,
    },
  });
}

describe("Compare review body", () => {
  beforeEach(() => {
    sessionStorage.clear();
    useAppState.setState({ activeDraft: null, username: undefined });
  });

  it("renders empty state under /experiments/:eid/compare with no active id", () => {
    renderAt("/experiments/7/compare");
    expect(screen.getByTestId("compare-page-body")).toBeInTheDocument();
    // No comparison selected → empty state placeholder
    expect(screen.getByTestId("compare-empty-state")).toBeInTheDocument();
  });

  it("reads :id from URL when present (review mode placeholder)", () => {
    renderAt("/experiments/7/compare/42");
    const page = screen.getByTestId("compare-page-body");
    expect(page).toBeInTheDocument();
    expect(page).toHaveAttribute("data-comparison-id", "42");
  });

  it("renders global listing scope under /compare/all", () => {
    renderAt("/compare/all");
    const body = screen.getByTestId("compare-page-body");
    expect(body).toBeInTheDocument();
    expect(body).toHaveAttribute("data-scope", "all");
  });
});

describe("Compare edit body", () => {
  beforeEach(() => {
    sessionStorage.clear();
    useAppState.setState({ activeDraft: null, username: undefined });
  });

  it("renders edit body under /experiments/:eid/compare/new", () => {
    renderAt("/experiments/7/compare/new");
    expect(screen.getByTestId("compare-page-edit")).toBeInTheDocument();
  });

  it("reads :id from URL under /experiments/:eid/compare/:id with an active draft", () => {
    seedDraftFor(42);
    renderAt("/experiments/7/compare/42");
    const edit = screen.getByTestId("compare-page-edit");
    expect(edit).toBeInTheDocument();
    expect(edit).toHaveAttribute("data-comparison-id", "42");
  });
});

// ── Regression: ResizeObserver re-attach when Skeleton swaps in (issue #51)

describe("Compare review mode — ResizeObserver", () => {
  let ROInstances: { observe: unknown; disconnect: unknown }[] = [];

  beforeEach(() => {
    sessionStorage.clear();
    useAppState.setState({ activeDraft: null, username: undefined });
    ROInstances = [];
    vi.stubGlobal("ResizeObserver", vi.fn((cb: ResizeObserverCallback) => {
      const inst = {
        observe: vi.fn((el: Element) => {
          // Give the element a non-zero clientHeight so the page-level
          // observer callback sets panelHeight > 0.
          Object.defineProperty(el, "clientHeight", { value: 400, configurable: true });
          cb([{ target: el, contentRect: { height: 400 } } as unknown as ResizeObserverEntry], inst as unknown as ResizeObserver);
        }),
        disconnect: vi.fn(),
      };
      ROInstances.push(inst);
      return inst;
    }));
  });

  it("attaches ResizeObserver after Skeleton lifts (plotLoading false)", async () => {
    const qc = new QueryClient({
      defaultOptions: { queries: { retry: false, gcTime: Infinity, staleTime: 0 }, mutations: { retry: false } },
    });
    // Do NOT seed the comparison query — let it start in isLoading so
    // Skeleton renders its fallback, then resolves asynchronously.
    // This exercises the bug: with deps: [] the effect runs while the
    // ref is null (Skeleton has no children) and never re-attaches.

    vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/comparisons/42") {
        return new Response(JSON.stringify({
          id: 42, title: "T", description: null, content_hash: "h",
          created_by: 1, created_at: "", updated_at: "",
          forked_from_id: null, forked_at_hash: null, forked_from_title: null,
          members: [{ id: 1, exposure_id: 100, display_order: 0, band_height: 1, y_offset: 0, normalization: "max", snapshot: null }],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/100/trace") {
        return new Response(JSON.stringify({ q: [0.1], I: [1], sigma: [0.01] }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });

    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Wait for the real plot column (not skeleton) to render.
    await waitFor(() => {
      expect(screen.queryByText("Loading comparison…")).toBeNull();
    });

    // The page-level ResizeObserver drives panelHeight state, which is
    // passed to MemberMetaGutter as style.height. Asserting on the gutter's
    // height pins to the page-level observer specifically (MultiTracePlot
    // has its own observer but does not affect panelHeight).
    const gutter = screen.getByTestId("member-meta-gutter");
    expect(gutter.style.height).not.toBe("0px");
  });
});

// ── Issues #61 + #52: cold-cache exposure/sample hydration ────────────────────
//
// Review mode loads `/experiments/:eid/compare/:id` directly. Before this fix
// the page only fetched comparison + traces — no exposure / sample rows. That
// left:
//   - `sampleIdFor` returning null for every member → ORPHAN_FALLBACK gray
//     for every trace stroke regardless of grouping mode (#61).
//   - `MemberMetaRow` falling back to `Exposure #N` instead of the
//     human-friendly `${sample.display_name || sample.name} · ${exposure.filename}` (#52).
// Both regressions resolve once review mode subscribes to each member's
// exposure row (and the matching sample row) so the cache populates.

describe("Compare review mode — cold-cache exposure + sample hydration (#61, #52)", () => {
  beforeEach(() => {
    sessionStorage.clear();
    useAppState.setState({ activeDraft: null, username: undefined });
    vi.stubGlobal("ResizeObserver", vi.fn(() => ({
      observe: vi.fn((el: Element) => {
        Object.defineProperty(el, "clientHeight", { value: 400, configurable: true });
      }),
      disconnect: vi.fn(),
    })));
  });

  function makeFetchMock() {
    return vi.fn(async (input: RequestInfo | URL) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/comparisons/42") {
        return new Response(JSON.stringify({
          id: 42, title: "T", description: null, content_hash: "h",
          created_by: 1, created_at: "2026-05-07T00:00:00",
          updated_at: "2026-05-07T00:00:00",
          forked_from_id: null, forked_at_hash: null, forked_from_title: null,
          members: [
            {
              id: 1, comparison_id: 42, exposure_id: 100, display_order: 0,
              band_height: 1, y_offset: 0, normalization: "none",
              color_override: null, label_override: null,
              q_window_min: null, q_window_max: null, peak_display: null,
              snapshot: null, is_stale: false, created_by: 1,
              created_at: "2026-05-07T00:00:00",
            },
            {
              id: 2, comparison_id: 42, exposure_id: 200, display_order: 1,
              band_height: 1, y_offset: 0, normalization: "none",
              color_override: null, label_override: null,
              q_window_min: null, q_window_max: null, peak_display: null,
              snapshot: null, is_stale: false, created_by: 1,
              created_at: "2026-05-07T00:00:00",
            },
          ],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/100") {
        return new Response(JSON.stringify({
          id: 100, sample_id: 11, filename: "run-A.dat", kind: "file",
          selected: true, status: "accepted", image_path: null,
          image_version: "", tags: [], sources: [],
          trace_hash: "h1", analysis_inputs_hash: "h1",
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/200") {
        return new Response(JSON.stringify({
          id: 200, sample_id: 22, filename: "run-B.dat", kind: "file",
          selected: true, status: "accepted", image_path: null,
          image_version: "", tags: [], sources: [],
          trace_hash: "h2", analysis_inputs_hash: "h2",
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/samples/11") {
        return new Response(JSON.stringify({
          id: 11, experiment_id: 7, display_name: "Lipid-A", name: "JC001",
          notes: null, tags: [],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/samples/22") {
        return new Response(JSON.stringify({
          id: 22, experiment_id: 7, display_name: "Lipid-B", name: "JC002",
          notes: null, tags: [],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/100/trace") {
        return new Response(JSON.stringify({ q: [0.1, 0.2], I: [1, 2], sigma: [0.01, 0.01] }), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      if (url === "/api/exposures/200/trace") {
        return new Response(JSON.stringify({ q: [0.1, 0.2], I: [3, 4], sigma: [0.01, 0.01] }), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response("not found", { status: 404 });
    }) as unknown as typeof fetch;
  }

  it("hydrates the per-member exposure cache so sampleIdFor resolves (#61)", async () => {
    vi.spyOn(global, "fetch").mockImplementation(makeFetchMock());
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
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Wait for the page to leave the cold-loading skeleton so all hooks
    // (comparison, traces, AND exposures) have had a chance to fire.
    await waitFor(() => {
      expect(screen.queryByText("Loading comparison…")).toBeNull();
    });

    // The bug: review mode never subscribed to per-member exposure rows, so
    // these cache reads stay `undefined` indefinitely, sampleIdFor returns
    // null, and every trace renders in ORPHAN_FALLBACK gray.
    await waitFor(() => {
      expect(qc.getQueryData(queryKeys.exposure(100))).toBeDefined();
      expect(qc.getQueryData(queryKeys.exposure(200))).toBeDefined();
    });
  });

  it("renders MemberMetaRow with sample label + filename, not 'Exposure #N' (#52)", async () => {
    vi.spyOn(global, "fetch").mockImplementation(makeFetchMock());
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
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    await waitFor(() => {
      expect(screen.queryByText("Loading comparison…")).toBeNull();
    });

    // After exposure + sample data lands, the gutter resolves human labels:
    // `${sample.display_name || sample.name} · ${exposure.filename}` (truncated).
    // The internal id form ("Exposure #100") is the regressed fallback we're moving away from.
    await waitFor(() => {
      const labels = screen.getAllByTestId("member-meta-label").map((n) => n.textContent ?? "");
      expect(labels.some((t) => t.includes("Lipid-A") && t.includes("run-A"))).toBe(true);
      expect(labels.some((t) => t.includes("Lipid-B") && t.includes("run-B"))).toBe(true);
      expect(labels.every((t) => !t.startsWith("Exposure #"))).toBe(true);
    });
  });
});
