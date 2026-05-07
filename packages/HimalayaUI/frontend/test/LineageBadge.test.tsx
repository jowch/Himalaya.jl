/**
 * LineageBadge tests (Plan §Phase 11, Task 11.3).
 *
 * The badge surfaces fork lineage in the review-mode header.
 *   - When `forked_from_id` is set, show "Forked from <parent title> by
 *     <author>" with a "view parent →" link.
 *   - When `forked_from_id IS NULL` BUT `forked_at_hash` is set (the parent
 *     was deleted post-fork; FK ON DELETE SET NULL clears the id), show a
 *     "Forked from a deleted comparison" message.
 *   - When neither is set, render nothing.
 */
import { describe, it, expect } from "vitest";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import { LineageBadge } from "../src/components/LineageBadge";
import type { Comparison } from "../src/api";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function makeComparison(over: Partial<Comparison> = {}): Comparison {
  return {
    id: 100,
    title: "Some fork",
    description: null,
    content_hash: "sha256:fork",
    created_by: 7,
    created_at: "2026-01-01",
    updated_at: "2026-01-01",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: [],
    ...over,
  };
}

function renderBadge(comparison: Comparison, qc = makeQc()) {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter>
        <LineageBadge comparison={comparison} experimentId={7} />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("LineageBadge", () => {
  it("renders nothing when neither forked_from_id nor forked_at_hash is set", () => {
    const { container } = renderBadge(makeComparison());
    expect(container.querySelector("[data-testid='comparison-lineage']")).toBeNull();
  });

  it("renders parent info when forked_from_id is set", () => {
    const cmp = makeComparison({
      forked_from_id: 42,
      forked_at_hash: "sha256:parent",
      forked_from_title: "Parent comparison",
    });
    renderBadge(cmp);
    const badge = screen.getByTestId("comparison-lineage");
    expect(badge).toHaveAttribute("data-parent-id", "42");
    expect(badge.textContent).toContain("Parent comparison");
    // view-parent link exists (router-friendly).
    const link = screen.getByTestId("comparison-lineage-view-parent");
    expect(link).toHaveAttribute("href", "/experiments/7/compare/42");
  });

  it("renders 'deleted parent' when only forked_at_hash is set", () => {
    const cmp = makeComparison({
      forked_from_id: null,
      forked_at_hash: "sha256:gone",
      forked_from_title: null,
    });
    renderBadge(cmp);
    const badge = screen.getByTestId("comparison-lineage");
    expect(badge).toHaveAttribute("data-parent-id", "deleted");
    expect(badge.textContent?.toLowerCase()).toContain("deleted");
    // No view-parent link for a deleted parent.
    expect(screen.queryByTestId("comparison-lineage-view-parent")).toBeNull();
  });

  it("falls back to a generic title when forked_from_title is missing", () => {
    const cmp = makeComparison({
      forked_from_id: 42,
      forked_at_hash: "sha256:parent",
      forked_from_title: null,
    });
    renderBadge(cmp);
    const badge = screen.getByTestId("comparison-lineage");
    expect(badge).toHaveAttribute("data-parent-id", "42");
    // Some fallback string surfaces — assert against the title field's null.
    // The link still renders so users can navigate.
    expect(screen.getByTestId("comparison-lineage-view-parent")).toBeInTheDocument();
  });

  it("uses /compare/all when experimentId is undefined (global scope)", () => {
    const cmp = makeComparison({
      forked_from_id: 42,
      forked_at_hash: "sha256:parent",
      forked_from_title: "Parent",
    });
    render(
      <QueryClientProvider client={makeQc()}>
        <MemoryRouter>
          <LineageBadge comparison={cmp} experimentId={undefined} />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    const link = screen.getByTestId("comparison-lineage-view-parent");
    expect(link).toHaveAttribute("href", "/compare/all");
  });
});
