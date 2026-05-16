/**
 * Fork affordance tests (Plan §Phase 11, Task 11.2; migrated by Compare UX
 * C-12).
 *
 * Compare UX C-12 replaced the review-mode header's inline `EditOrForkButton`
 * with the `CompareToolbar` overflow menu. The Edit gesture is deferred to the
 * SavePill in C-12's successor (C-14); only the Fork action survives in the
 * toolbar, and it is no longer author-gated — anyone can fork from the
 * overflow menu.
 *
 * This file verifies the Fork action wiring through the new toolbar:
 *   - Fork lives in the `⋯ more` overflow menu and is always present.
 *   - Fork click creates a brand-new draft (`id === undefined`) carrying
 *     `forkedFromId` + `forkedAtHash` from the parent, and navigates to
 *     /experiments/:eid/compare/new.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Compare } from "../src/pages/Compare";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";

// MultiTracePlot pulls Observable Plot — stub it to avoid jsdom weirdness.
vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    (el as unknown as { scale: (n: string) => unknown }).scale = (n) =>
      n === "x"
        ? { invert: (px: number) => px / 100, apply: (q: number) => q * 100 }
        : undefined;
    return el;
  }),
  line: vi.fn(() => ({})),
  dot:  vi.fn(() => ({})),
  text: vi.fn(() => ({})),
  link: vi.fn(() => ({})),
}));

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
    id: 42,
    title: "T",
    description: "Some description",
    content_hash: "sha256:parent",
    created_by: 7,
    created_at: "2026-01-01",
    updated_at: "2026-01-01",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: [
      {
        id: 100, comparison_id: 42, exposure_id: null,
        display_order: 0, band_height: 1, y_offset: 0,
        normalization: "qwindow",
        color_override: null, label_override: null,
        q_window_min: null, q_window_max: null, peak_display: null,
        snapshot: {
          effective_peaks: [],
          confirmed_index: null,
          analysis_inputs_hash: "h",
        },
        is_stale: false, created_by: null, created_at: null,
      },
    ],
    ...over,
  };
}

function LocationSpy() {
  const loc = useLocation();
  return <div data-testid="current-location">{loc.pathname}</div>;
}

function renderReview(qc: QueryClient) {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
        <LocationSpy />
        <Routes>
          <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

/** Open the toolbar overflow menu and return its `Fork` button. */
async function openForkMenuItem(): Promise<HTMLElement> {
  const more = await waitFor(() => screen.getByTestId("compare-toolbar-more"));
  fireEvent.click(more);
  return waitFor(() => screen.getByRole("button", { name: "Fork" }));
}

describe("Fork affordance — toolbar overflow menu", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
    useAppState.getState().discardDraft();
  });

  it("Fork is reachable from the toolbar overflow menu", async () => {
    const qc = makeQc();
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));

    renderReview(qc);
    const forkItem = await openForkMenuItem();
    expect(forkItem).toBeInTheDocument();
  });

  it("Fork is available regardless of authorship (non-author)", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 99, username: "bob", first_name: null, last_name: null },
    ]);
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));
    useAppState.setState({ username: "bob" });

    renderReview(qc);
    const forkItem = await openForkMenuItem();
    expect(forkItem).toBeInTheDocument();
  });
});

describe("Fork affordance — action", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
    useAppState.getState().discardDraft();
  });

  it("Fork click creates a fresh draft with parent lineage and navigates to /new", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 99, username: "bob", first_name: null, last_name: null },
    ]);
    qc.setQueryData(["comparison", 42] as const, makeComparison({
      created_by: 7,
      title: "Parent comparison",
      description: "Parent description",
    }));
    useAppState.setState({ username: "bob" });

    renderReview(qc);
    const forkItem = await openForkMenuItem();
    fireEvent.click(forkItem);
    await waitFor(() =>
      expect(screen.getByTestId("current-location").textContent)
        .toBe("/experiments/7/compare/new"),
    );
    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    // Brand-new draft (no server id; no baseHash).
    expect(draft?.id).toBeUndefined();
    expect(draft?.baseHash).toBeUndefined();
    // Lineage rides through.
    expect(draft?.forkedFromId).toBe(42);
    expect(draft?.forkedAtHash).toBe("sha256:parent");
    // Title prefilled from parent.
    expect(draft?.title).toBe("Fork of Parent comparison");
    // Members projected from parent — but server id dropped (each member
    // becomes a fresh INSERT under the new comparison).
    expect(draft?.members.length).toBe(1);
    expect(draft?.members[0]?.id).toBeUndefined();
  });
});
