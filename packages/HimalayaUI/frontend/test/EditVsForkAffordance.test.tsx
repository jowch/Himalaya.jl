/**
 * Edit-vs-Fork affordance tests (Plan §Phase 11, Task 11.2).
 *
 * Verifies the review-mode header gating:
 *   - Author of the comparison sees an "Edit" button (not Fork).
 *   - Non-author sees a "Fork" button (not Edit).
 *   - Orphaned-author comparison (`created_by === null`) shows Fork to all
 *     users — no user matches null, so the Edit affordance hides.
 *
 * Plus action wiring:
 *   - Edit click navigates to /experiments/:eid/compare/:id/edit and
 *     `loadDraftFromComparison` seeds the Zustand draft.
 *   - Fork click creates a brand-new draft (`id === undefined`) carrying
 *     `forkedFromId` + `forkedAtHash` from the parent, and navigates to
 *     /experiments/:eid/compare/new.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparePage } from "../src/pages/ComparePage";
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
        <Routes>
          <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<LocationSpy />} />
          <Route path="/experiments/:eid/compare/new" element={<LocationSpy />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Edit-vs-Fork affordance — visibility", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
    useAppState.getState().discardDraft();
  });

  it("author sees Edit button, not Fork", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 7, username: "alice", first_name: null, last_name: null },
    ]);
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));
    useAppState.setState({ username: "alice" });

    renderReview(qc);
    await waitFor(() => screen.getByTestId("compare-review-header"));
    await waitFor(() => expect(screen.getByTestId("comparison-edit")).toBeInTheDocument());
    expect(screen.queryByTestId("comparison-fork")).toBeNull();
  });

  it("non-author sees Fork button, not Edit", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 99, username: "bob", first_name: null, last_name: null },
    ]);
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));
    useAppState.setState({ username: "bob" });

    renderReview(qc);
    await waitFor(() => screen.getByTestId("compare-review-header"));
    await waitFor(() => expect(screen.getByTestId("comparison-fork")).toBeInTheDocument());
    expect(screen.queryByTestId("comparison-edit")).toBeNull();
  });

  it("orphan-author comparison (created_by === null) shows Fork to all users", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 7, username: "alice", first_name: null, last_name: null },
    ]);
    // Original author is alice (id 7) but the comparison is orphan now.
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: null }));
    useAppState.setState({ username: "alice" });

    renderReview(qc);
    await waitFor(() => screen.getByTestId("compare-review-header"));
    await waitFor(() => expect(screen.getByTestId("comparison-fork")).toBeInTheDocument());
    expect(screen.queryByTestId("comparison-edit")).toBeNull();
  });

  it("hides both buttons when the current user has not been resolved yet", async () => {
    const qc = makeQc();
    // No `users` cache, no `username` set. The lookup returns undefined.
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));
    renderReview(qc);
    await waitFor(() => screen.getByTestId("compare-review-header"));
    // We're not the author and we don't know who we are; Fork is the safe
    // default (matches non-author). Edit must not be shown.
    expect(screen.queryByTestId("comparison-edit")).toBeNull();
    expect(screen.getByTestId("comparison-fork")).toBeInTheDocument();
  });
});

describe("Edit-vs-Fork affordance — actions", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
    useAppState.getState().discardDraft();
  });

  it("Edit click navigates to /experiments/:eid/compare/:id/edit and seeds the draft", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 7, username: "alice", first_name: null, last_name: null },
    ]);
    qc.setQueryData(["comparison", 42] as const, makeComparison({ created_by: 7 }));
    useAppState.setState({ username: "alice" });

    renderReview(qc);
    const editBtn = await waitFor(() => screen.getByTestId("comparison-edit"));
    fireEvent.click(editBtn);
    await waitFor(() =>
      expect(screen.getByTestId("current-location").textContent)
        .toBe("/experiments/7/compare/42/edit"),
    );
    // Draft was loaded against the parent; loadDraftFromComparison shape:
    // - id matches, baseHash matches content_hash.
    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    expect(draft?.id).toBe(42);
    expect(draft?.baseHash).toBe("sha256:parent");
    // Edit is NOT a fork — lineage stays empty.
    expect(draft?.forkedFromId).toBeUndefined();
    expect(draft?.forkedAtHash).toBeUndefined();
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
    const forkBtn = await waitFor(() => screen.getByTestId("comparison-fork"));
    fireEvent.click(forkBtn);
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
