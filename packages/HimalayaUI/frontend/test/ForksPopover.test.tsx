/**
 * ForksPopover tests (Plan §Phase 11, Task 11.3).
 *
 * - "Forks (N) →" trigger surfaces the count from `useComparisonForks`.
 * - Click toggles a popover listing the forks (one row per fork, each
 *   linking to the fork's review page).
 * - Empty state: shows "No forks yet" when the fork list is empty.
 */
import { describe, it, expect } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ForksPopover } from "../src/components/ForksPopover";
import type { ComparisonSummary } from "../src/api";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function makeFork(over: Partial<ComparisonSummary> = {}): ComparisonSummary {
  return {
    id: 200,
    title: "Some fork",
    description: null,
    content_hash: "sha256:fork",
    created_by: 7,
    created_at: "2026-01-02",
    updated_at: "2026-01-02",
    forked_from_id: 100,
    forked_at_hash: "sha256:parent",
    ...over,
  };
}

function renderPopover(opts: {
  parentId: number;
  experimentId?: number | undefined;
  forks: ComparisonSummary[];
}) {
  const qc = makeQc();
  qc.setQueryData(["comparison", opts.parentId, "forks"] as const, opts.forks);
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter>
        <ForksPopover comparisonId={opts.parentId} experimentId={opts.experimentId} />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ForksPopover", () => {
  it("trigger renders 'Forks (N)' with the cached count", () => {
    renderPopover({
      parentId: 100,
      experimentId: 7,
      forks: [makeFork({ id: 200 }), makeFork({ id: 201 })],
    });
    const trigger = screen.getByTestId("comparison-forks-trigger");
    expect(trigger).toBeInTheDocument();
    expect(trigger.textContent).toMatch(/Forks \(2\)/);
  });

  it("trigger renders zero count when there are no forks", () => {
    renderPopover({ parentId: 100, experimentId: 7, forks: [] });
    const trigger = screen.getByTestId("comparison-forks-trigger");
    expect(trigger.textContent).toMatch(/Forks \(0\)/);
  });

  it("popover opens on trigger click and lists fork rows", () => {
    renderPopover({
      parentId: 100,
      experimentId: 7,
      forks: [
        makeFork({ id: 200, title: "First fork" }),
        makeFork({ id: 201, title: "Second fork" }),
      ],
    });
    // Closed initially.
    expect(screen.queryByTestId("comparison-forks-popover")).toBeNull();

    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    const popover = screen.getByTestId("comparison-forks-popover");
    expect(popover).toBeInTheDocument();
    expect(screen.getAllByTestId("comparison-forks-row")).toHaveLength(2);
    expect(popover.textContent).toContain("First fork");
    expect(popover.textContent).toContain("Second fork");
  });

  it("each fork row links to its review page in the right scope", () => {
    renderPopover({
      parentId: 100,
      experimentId: 7,
      forks: [makeFork({ id: 200, title: "Linked fork" })],
    });
    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    const row = screen.getByTestId("comparison-forks-row");
    const link = row.querySelector("a");
    expect(link).not.toBeNull();
    expect(link).toHaveAttribute("href", "/experiments/7/compare/200");
  });

  it("renders 'No forks yet' empty state when the list is empty", () => {
    renderPopover({ parentId: 100, experimentId: 7, forks: [] });
    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    const popover = screen.getByTestId("comparison-forks-popover");
    expect(popover.textContent?.toLowerCase()).toContain("no forks yet");
  });

  it("global scope deep-links forks to /compare/all/:forkId", () => {
    renderPopover({
      parentId: 100,
      experimentId: undefined,
      forks: [makeFork({ id: 200, title: "Fork" })],
    });
    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    const link = screen.getByTestId("comparison-forks-row").querySelector("a");
    // Phase 4 follow-up: the global scope has /compare/all/:id, so each
    // fork row deep-links directly rather than punting the user to the
    // empty global listing.
    expect(link).toHaveAttribute("href", "/compare/all/200");
  });

  it("clicking the trigger again closes the popover", () => {
    renderPopover({ parentId: 100, experimentId: 7, forks: [] });
    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    expect(screen.getByTestId("comparison-forks-popover")).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
    expect(screen.queryByTestId("comparison-forks-popover")).toBeNull();
  });
});
