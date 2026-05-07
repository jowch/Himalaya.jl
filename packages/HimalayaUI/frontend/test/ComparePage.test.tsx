/**
 * Page shells (Plan §Phase 4, Task 4.1).
 *
 * Verifies that ComparePage and ComparePageEdit render under the expected
 * routes and read URL params correctly. Behaviour exercised by later tasks
 * (sidebar, draft state, save flow) lives in their own test files.
 */
import { describe, it, expect } from "vitest";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
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
