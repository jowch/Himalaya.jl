/**
 * ComparisonSidebar (Plan §Phase 4, Task 4.2).
 *
 * Asserts:
 * - Sort order (most-recent first by `updated_at`).
 * - Scope toggle: clicking "All experiments" navigates to /compare/all.
 * - Search filter narrows the rendered list.
 * - Empty state when no comparisons.
 * - Active comparison highlighted via `data-active="true"`.
 *
 * Tests use a fresh QueryClient + MemoryRouter, seed the comparison list
 * via `qc.setQueryData`, and inspect rendered DOM for stable selectors.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, within, waitFor, fireEvent } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { useAppState } from "../src/state";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparisonSidebar } from "../src/components/ComparisonSidebar";
import { queryKeys } from "../src/queries";
import type { ComparisonSummary } from "../src/api";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderSidebar(opts: {
  qc: QueryClient;
  scope: "experiment" | "all";
  experimentId: number | undefined;
  activeId?: number | undefined;
  initialPath?: string;
  /** Probe component renders the URL pathname for assertions. */
  probe?: boolean;
}) {
  const path = opts.initialPath
    ?? (opts.scope === "experiment"
      ? `/experiments/${opts.experimentId}/compare`
      : "/compare/all");
  return render(
    <QueryClientProvider client={opts.qc}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route
            path="/experiments/:eid/compare"
            element={
              <>
                <ComparisonSidebar
                  experimentId={opts.experimentId}
                  scope={opts.scope}
                  activeComparisonId={opts.activeId}
                />
                {opts.probe && <PathProbe />}
              </>
            }
          />
          <Route
            path="/experiments/:eid/compare/:id"
            element={
              <>
                <ComparisonSidebar
                  experimentId={opts.experimentId}
                  scope={opts.scope}
                  activeComparisonId={opts.activeId}
                />
                {opts.probe && <PathProbe />}
              </>
            }
          />
          <Route
            path="/compare/all"
            element={
              <>
                <ComparisonSidebar
                  experimentId={undefined}
                  scope="all"
                  activeComparisonId={opts.activeId}
                />
                {opts.probe && <PathProbe />}
              </>
            }
          />
          <Route
            path="/experiments/:eid/compare/new"
            element={
              <>
                <ComparisonSidebar
                  experimentId={opts.experimentId}
                  scope={opts.scope}
                  activeComparisonId={opts.activeId}
                />
                {opts.probe && <PathProbe />}
              </>
            }
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

import { useLocation } from "react-router-dom";
function PathProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="path-probe">{loc.pathname}</div>;
}

const ROW_OLD: ComparisonSummary = {
  id: 1,
  title: "Old comparison",
  description: null,
  content_hash: "h1",
  created_by: 1,
  created_at: "2026-04-01T00:00:00Z",
  updated_at: "2026-04-01T00:00:00Z",
  forked_from_id: null,
  forked_at_hash: null,
};

const ROW_MID: ComparisonSummary = {
  ...ROW_OLD,
  id: 2,
  title: "Mid comparison",
  updated_at: "2026-04-15T00:00:00Z",
};

const ROW_NEW: ComparisonSummary = {
  ...ROW_OLD,
  id: 3,
  title: "Recent comparison",
  updated_at: "2026-05-01T00:00:00Z",
};

describe("ComparisonSidebar", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = makeQc();
  });

  it("renders empty state when no comparisons", () => {
    qc.setQueryData(queryKeys.comparisons(7), []);
    renderSidebar({ qc, scope: "experiment", experimentId: 7 });
    expect(screen.getByTestId("comparison-sidebar-empty")).toBeInTheDocument();
  });

  it("sorts rows most-recent first by updated_at", () => {
    qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD, ROW_NEW, ROW_MID]);
    renderSidebar({ qc, scope: "experiment", experimentId: 7 });
    const rows = screen.getAllByTestId("comparison-list-item");
    expect(rows).toHaveLength(3);
    expect(rows[0]).toHaveAttribute("data-comparison-id", "3");
    expect(rows[1]).toHaveAttribute("data-comparison-id", "2");
    expect(rows[2]).toHaveAttribute("data-comparison-id", "1");
  });

  it("highlights the active comparison", () => {
    qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD, ROW_NEW]);
    renderSidebar({ qc, scope: "experiment", experimentId: 7, activeId: 3 });
    const rows = screen.getAllByTestId("comparison-list-item");
    const active = rows.find((r) => r.getAttribute("data-comparison-id") === "3")!;
    const inactive = rows.find((r) => r.getAttribute("data-comparison-id") === "1")!;
    expect(active).toHaveAttribute("data-active", "true");
    expect(inactive).not.toHaveAttribute("data-active", "true");
  });

  it("scope toggle exposes both options with `this` selected by default", () => {
    qc.setQueryData(queryKeys.comparisons(7), []);
    renderSidebar({ qc, scope: "experiment", experimentId: 7 });
    const toggle = screen.getByTestId("comparison-scope-toggle");
    expect(toggle).toHaveAttribute("data-scope", "this");
    expect(within(toggle).getByText(/this experiment/i)).toBeInTheDocument();
    expect(within(toggle).getByText(/all experiments/i)).toBeInTheDocument();
  });

  it("clicking 'All experiments' navigates to /compare/all", async () => {
    const user = userEvent.setup();
    qc.setQueryData(queryKeys.comparisons(7), []);
    qc.setQueryData(queryKeys.comparisons("all"), []);
    renderSidebar({
      qc, scope: "experiment", experimentId: 7, probe: true,
    });
    await user.click(screen.getByTestId("comparison-scope-all"));
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/compare/all");
  });

  it("search box filters the rendered list", async () => {
    const user = userEvent.setup();
    qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD, ROW_MID, ROW_NEW]);
    renderSidebar({ qc, scope: "experiment", experimentId: 7 });
    const search = screen.getByTestId("comparison-sidebar-search");
    await user.type(search, "Recent");
    const rows = screen.getAllByTestId("comparison-list-item");
    expect(rows).toHaveLength(1);
    expect(rows[0]).toHaveAttribute("data-comparison-id", "3");
  });

  it("'+ New' button navigates to /experiments/:eid/compare/new", async () => {
    const user = userEvent.setup();
    qc.setQueryData(queryKeys.comparisons(7), []);
    renderSidebar({
      qc, scope: "experiment", experimentId: 7, probe: true,
    });
    await user.click(screen.getByTestId("comparison-new"));
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/new");
  });

  // ─── Pin tests (Phase 13, Task 13.2) ──────────────────────────────────────

  describe("pins", () => {
    beforeEach(() => {
      // The pin hooks read `useAppState((s) => s.username)`. Set it so the
      // GET /api/users/me/comparison-pins query is `enabled`.
      useAppState.setState({ username: "alice" });
    });

    it("pinned items render at top of list, marked with data-pinned", () => {
      qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD, ROW_MID, ROW_NEW]);
      // Pin the OLDEST row → it should float to the top regardless of
      // updated_at (pins beat updated_at).
      qc.setQueryData(queryKeys.comparisonPins, [1]);
      renderSidebar({ qc, scope: "experiment", experimentId: 7 });
      const rows = screen.getAllByTestId("comparison-list-item");
      expect(rows[0]).toHaveAttribute("data-comparison-id", "1");
      expect(rows[0]).toHaveAttribute("data-pinned", "true");
      // The other two stay in updated_at order.
      expect(rows[1]).toHaveAttribute("data-comparison-id", "3");
      expect(rows[2]).toHaveAttribute("data-comparison-id", "2");
    });

    it("multiple pins preserve API order (most-recently-pinned first)", () => {
      qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD, ROW_MID, ROW_NEW]);
      // API returns pins in most-recent-pin-first; the sidebar must
      // preserve that order (NOT re-sort by updated_at).
      qc.setQueryData(queryKeys.comparisonPins, [2, 3]); // 2 pinned more recently
      renderSidebar({ qc, scope: "experiment", experimentId: 7 });
      const rows = screen.getAllByTestId("comparison-list-item");
      expect(rows[0]).toHaveAttribute("data-comparison-id", "2");
      expect(rows[1]).toHaveAttribute("data-comparison-id", "3");
      expect(rows[2]).toHaveAttribute("data-comparison-id", "1");
    });

    it("clicking the pin toggle on an unpinned row POSTs /pin", async () => {
      qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD]);
      qc.setQueryData(queryKeys.comparisonPins, []);
      const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
        const url = typeof input === "string" ? input : String(input);
        if (url === "/api/comparisons/1/pin" && init?.method === "POST") {
          return new Response(
            JSON.stringify({ comparison_id: 1, pinned: true }),
            { status: 200, headers: { "Content-Type": "application/json" } },
          );
        }
        return new Response("[]", {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      });
      renderSidebar({ qc, scope: "experiment", experimentId: 7 });
      const toggle = await screen.findByTestId("comparison-pin-toggle");
      expect(toggle).toHaveAttribute("data-pinned", "false");
      // fireEvent over userEvent because the actionability checks userEvent
      // performs (opacity, pointer-events, etc.) trip on the hover-revealed
      // pin toggle's `opacity-0` Tailwind class even though JSDOM doesn't
      // actually compute styles. fireEvent dispatches the synthetic event
      // directly which is the right contract for testing handler invocation.
      fireEvent.click(toggle);
      await waitFor(() => {
        const calls = fetchSpy.mock.calls.map((c) =>
          [typeof c[0] === "string" ? c[0] : String(c[0]), (c[1] as RequestInit | undefined)?.method] as const,
        );
        expect(calls).toContainEqual(["/api/comparisons/1/pin", "POST"]);
      });
    });

    it("clicking the pin toggle on a pinned row DELETEs /pin", async () => {
      qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD]);
      qc.setQueryData(queryKeys.comparisonPins, [1]);
      const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
        const url = typeof input === "string" ? input : String(input);
        if (url === "/api/comparisons/1/pin" && init?.method === "DELETE") {
          return new Response(
            JSON.stringify({ comparison_id: 1, pinned: false }),
            { status: 200, headers: { "Content-Type": "application/json" } },
          );
        }
        return new Response("[]", {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      });
      renderSidebar({ qc, scope: "experiment", experimentId: 7 });
      const toggle = await screen.findByTestId("comparison-pin-toggle");
      expect(toggle).toHaveAttribute("data-pinned", "true");
      fireEvent.click(toggle);
      await waitFor(() => {
        const calls = fetchSpy.mock.calls.map((c) =>
          [typeof c[0] === "string" ? c[0] : String(c[0]), (c[1] as RequestInit | undefined)?.method] as const,
        );
        expect(calls).toContainEqual(["/api/comparisons/1/pin", "DELETE"]);
      });
    });

    it("pin toggle click does not bubble to the row navigation", () => {
      qc.setQueryData(queryKeys.comparisons(7), [ROW_OLD]);
      qc.setQueryData(queryKeys.comparisonPins, []);
      vi.spyOn(global, "fetch").mockResolvedValue(
        new Response("{}", { status: 200, headers: { "Content-Type": "application/json" } }),
      );
      renderSidebar({
        qc, scope: "experiment", experimentId: 7, probe: true,
      });
      const initialPath = screen.getByTestId("path-probe").textContent;
      fireEvent.click(screen.getByTestId("comparison-pin-toggle"));
      // Path must remain unchanged — clicking the pin button should not
      // also trigger row navigation via event bubbling.
      expect(screen.getByTestId("path-probe")).toHaveTextContent(initialPath ?? "");
    });
  });
});
