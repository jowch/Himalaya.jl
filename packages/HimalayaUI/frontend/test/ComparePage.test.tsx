/**
 * Page shells (Plan §Phase 4, Task 4.1).
 *
 * Verifies that ComparePage and ComparePageEdit render under the expected
 * routes and read URL params correctly. Behaviour exercised by later tasks
 * (sidebar, draft state, save flow) lives in their own test files.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor } from "@testing-library/react";
import { ComparePage } from "../src/pages/ComparePage";
import { ComparePageEdit } from "../src/pages/ComparePageEdit";

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
          <Route path="/experiments/:eid/compare" element={<ComparePage />} />
          <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          <Route path="/experiments/:eid/compare/new" element={<ComparePageEdit />} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<ComparePageEdit />} />
          <Route path="/compare/all" element={<ComparePage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ComparePage shell", () => {
  it("renders empty state under /experiments/:eid/compare with no active id", () => {
    renderAt("/experiments/7/compare");
    expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    // No comparison selected → empty state placeholder
    expect(screen.getByTestId("compare-empty-state")).toBeInTheDocument();
  });

  it("reads :id from URL when present (review mode placeholder)", () => {
    renderAt("/experiments/7/compare/42");
    const page = screen.getByTestId("compare-page");
    expect(page).toBeInTheDocument();
    expect(page).toHaveAttribute("data-comparison-id", "42");
  });

  it("renders global listing scope under /compare/all", () => {
    renderAt("/compare/all");
    expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    expect(screen.getByTestId("compare-page")).toHaveAttribute("data-scope", "all");
  });
});

describe("ComparePageEdit shell", () => {
  it("renders edit shell under /experiments/:eid/compare/new", () => {
    renderAt("/experiments/7/compare/new");
    expect(screen.getByTestId("compare-page-edit")).toBeInTheDocument();
  });

  it("reads :id from URL under /experiments/:eid/compare/:id/edit", () => {
    renderAt("/experiments/7/compare/42/edit");
    const edit = screen.getByTestId("compare-page-edit");
    expect(edit).toBeInTheDocument();
    expect(edit).toHaveAttribute("data-comparison-id", "42");
  });
});

// ── Regression: ResizeObserver re-attach when Skeleton swaps in (issue #51)

describe("ComparePage review mode — ResizeObserver", () => {
  let ROInstances: { observe: unknown; disconnect: unknown }[] = [];

  beforeEach(() => {
    ROInstances = [];
    vi.stubGlobal("ResizeObserver", vi.fn((cb: ResizeObserverCallback) => {
      const inst = {
        observe: vi.fn((el: Element) => { cb([{ target: el, contentRect: { height: 400 } } as unknown as ResizeObserverEntry], inst as unknown as ResizeObserver); }),
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
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Wait for the real plot column (not skeleton) to render.
    await waitFor(() => {
      expect(screen.queryByText("Loading comparison…")).toBeNull();
    });

    // With the bug (deps: []) the effect runs once while the ref is null
    // and the observer is never created. With the fix (deps: [plotLoading])
    // the effect re-runs when Skeleton swaps in children and the observer
    // attaches to the real DOM element.
    expect(ROInstances.length).toBeGreaterThanOrEqual(1);
    expect(ROInstances[0]!.observe).toHaveBeenCalled();
  });
});
