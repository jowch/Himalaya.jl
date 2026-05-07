/**
 * ConflictModal — Compare page 409 conflict resolution (Plan §Phase 12).
 *
 * Asserts:
 * - Hidden when `pendingConflict === null`.
 * - Side-by-side panels render the right titles + member counts.
 * - Discard: clears draft + clears slot + navigates to server's review.
 * - Overwrite: re-submits with `expected_content_hash` = server's hash;
 *   on success clears modal + navigates.
 * - Second-409 race: a re-submit returning a NEWER `current_state` keeps
 *   the modal open and updates the rendered server-state panel.
 * - Fork: starts a fork-flavored draft from the server's current state +
 *   navigates to the new-draft route.
 * - Esc closes the modal without committing any action (draft preserved).
 * - Outside-click closes without acting.
 * - Focus trap: Tab cycles within the dialog.
 *
 * The test wires the modal inside a MemoryRouter at
 * /experiments/7/compare/42/edit so :eid extraction can resolve to 7.
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { render, screen, act, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { ConflictModal } from "../src/components/ConflictModal";
import { useAppState, LS_KEY } from "../src/state";
import { ConflictError } from "../src/api";
import type { Comparison } from "../src/api";
import { COMPARE_DRAFT_KEY } from "../src/lib/comparison/draft";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function buildComparison(opts: {
  id?: number;
  hash?: string;
  title?: string;
  members?: number;
  description?: string | null;
  updated_at?: string | null;
}): Comparison {
  const {
    id = 42, hash = "sha256:server", title = "Server title",
    members = 2, description = null, updated_at = "2026-05-06T12:00:00Z",
  } = opts;
  return {
    id, title,
    description,
    content_hash: hash,
    created_by: 1,
    created_at: "2026-05-01T00:00:00Z",
    updated_at,
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: Array.from({ length: members }, (_, i) => ({
      id: 100 + i,
      comparison_id: id,
      exposure_id: 200 + i,
      display_order: i,
      band_height: 1,
      y_offset: 0,
      normalization: "none",
      color_override: null,
      label_override: null,
      q_window_min: null,
      q_window_max: null,
      peak_display: null,
      snapshot: {
        effective_peaks: [],
        confirmed_index: null,
        analysis_inputs_hash: "",
      },
      is_stale: false,
      created_by: 1,
      created_at: "2026-05-01T00:00:00Z",
    })),
  };
}

function PathProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="path-probe">{loc.pathname}</div>;
}

function renderModal(opts?: {
  initialPath?: string;
}): { qc: QueryClient } {
  const qc = makeQc();
  const initialPath = opts?.initialPath ?? "/experiments/7/compare/42/edit";
  const Wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <Routes>
          <Route path="*" element={
            <>
              {children}
              <PathProbe />
            </>
          } />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>
  );
  render(<ConflictModal />, { wrapper: Wrapper });
  return { qc };
}

function seedDraft(opts: { title?: string; members?: number; id?: number; baseHash?: string }) {
  const { title = "Local draft title", members = 3, id = 42, baseHash = "sha256:local" } = opts;
  useAppState.setState({
    activeDraft: {
      id,
      baseHash,
      title,
      description: "",
      members: Array.from({ length: members }, (_, i) => ({
        id: undefined,
        exposure_id: 300 + i,
        display_order: i,
        band_height: 1,
        y_offset: 0,
        normalization: "none",
        color_override: undefined,
        label_override: undefined,
        q_window_min: undefined,
        q_window_max: undefined,
        peak_display: undefined,
        snapshot: {
          effective_peaks: [],
          confirmed_index: null,
          analysis_inputs_hash: "",
        },
      })),
      forkedFromId: undefined,
      forkedAtHash: undefined,
    },
  });
}

beforeEach(() => {
  localStorage.clear();
  sessionStorage.clear();
  useAppState.setState({
    pendingConflict: null,
    activeDraft: null,
    username: "alice",
  });
});

afterEach(() => {
  vi.restoreAllMocks();
  localStorage.removeItem(LS_KEY);
  sessionStorage.removeItem(COMPARE_DRAFT_KEY);
});

describe("<ConflictModal>", () => {
  it("renders nothing when pendingConflict is null", () => {
    renderModal();
    expect(screen.queryByTestId("conflict-modal")).not.toBeInTheDocument();
  });

  it("renders side-by-side panels when conflict is set", () => {
    const server = buildComparison({ id: 42, title: "Server title", members: 2 });
    seedDraft({ title: "Local title", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    expect(screen.getByTestId("conflict-modal")).toBeInTheDocument();
    const serverPanel = screen.getByTestId("conflict-panel-server");
    const localPanel  = screen.getByTestId("conflict-panel-local");
    expect(within(serverPanel).getByTestId("conflict-panel-server-title")).toHaveTextContent("Server title");
    expect(within(serverPanel).getByTestId("conflict-panel-server-members")).toHaveTextContent("2");
    expect(within(localPanel).getByTestId("conflict-panel-local-title")).toHaveTextContent("Local title");
    expect(within(localPanel).getByTestId("conflict-panel-local-members")).toHaveTextContent("5");
  });

  it("Discard clears draft + clears slot + navigates to server review", async () => {
    const user = userEvent.setup();
    const server = buildComparison({ id: 42 });
    seedDraft({ title: "Local title", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    await user.click(screen.getByTestId("conflict-discard"));
    expect(useAppState.getState().pendingConflict).toBeNull();
    expect(useAppState.getState().activeDraft).toBeNull();
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42");
    expect(screen.queryByTestId("conflict-modal")).not.toBeInTheDocument();
  });

  it("Overwrite re-submits with expected_content_hash = server's current_hash, then closes + navigates on success", async () => {
    const user = userEvent.setup();
    const server = buildComparison({ id: 42, hash: "sha256:server-v1" });
    seedDraft({ title: "Local title", members: 1, id: 42, baseHash: "sha256:stale" });

    const captured: { url: string; body: unknown }[] = [];
    globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = typeof input === "string" ? input : String(input);
      const body = typeof init?.body === "string" ? JSON.parse(init.body) : undefined;
      captured.push({ url, body });
      return new Response(JSON.stringify(buildComparison({
        id: 42, hash: "sha256:after-overwrite", title: "Local title",
      })), { status: 200, headers: { "Content-Type": "application/json" } });
    }) as typeof fetch;

    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server-v1", server),
      );
    });
    await user.click(screen.getByTestId("conflict-overwrite"));

    await waitFor(() => expect(captured).toHaveLength(1));
    expect(captured[0]!.url).toBe("/api/comparisons/42/submit");
    expect(captured[0]!.body).toMatchObject({
      title: "Local title",
      expected_content_hash: "sha256:server-v1",
    });

    // Modal closes and navigates after success
    await waitFor(() => {
      expect(useAppState.getState().pendingConflict).toBeNull();
    });
    await waitFor(() => {
      expect(screen.queryByTestId("conflict-modal")).not.toBeInTheDocument();
    });
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42");
  });

  it("Overwrite on second-409 race: modal stays open with NEW server state", async () => {
    const user = userEvent.setup();
    const server1 = buildComparison({
      id: 42, hash: "sha256:v1", title: "Server v1", members: 2,
    });
    seedDraft({ title: "Local title", members: 1, id: 42, baseHash: "sha256:base" });

    // Second-409: the re-submit hits the server, but the server's hash has
    // moved on AGAIN between the first 409 and the re-submit.
    const server2 = buildComparison({
      id: 42, hash: "sha256:v2", title: "Server v2", members: 4,
    });
    globalThis.fetch = (async () => {
      return new Response(JSON.stringify({
        error: "conflict",
        current_hash: "sha256:v2",
        current_state: server2,
      }), { status: 409, headers: { "Content-Type": "application/json" } });
    }) as typeof fetch;

    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:v1", server1),
      );
    });
    expect(screen.getByTestId("conflict-panel-server-title")).toHaveTextContent("Server v1");

    await user.click(screen.getByTestId("conflict-overwrite"));

    // Wait for the second 409 to land + bridge into pendingConflict.
    await waitFor(() => {
      const slot = useAppState.getState().pendingConflict;
      expect(slot?.current_hash).toBe("sha256:v2");
    });

    // Modal still mounted; server-state panel reflects the NEW state.
    expect(screen.getByTestId("conflict-modal")).toBeInTheDocument();
    expect(screen.getByTestId("conflict-panel-server-title")).toHaveTextContent("Server v2");
    expect(screen.getByTestId("conflict-panel-server-members")).toHaveTextContent("4");
    // Path unchanged (no navigation on second 409)
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42/edit");
  });

  it("Fork starts a fork-draft from server's current state and navigates to new-draft", async () => {
    const user = userEvent.setup();
    const server = buildComparison({
      id: 42, hash: "sha256:server", title: "Server title", members: 2,
    });
    seedDraft({ title: "Local title", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    await user.click(screen.getByTestId("conflict-fork"));

    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    // fromComparisonAsFork sets forkedFromId/Hash from the SERVER state, and
    // member ids are dropped (each becomes an INSERT under the new comparison).
    expect(draft!.forkedFromId).toBe(42);
    expect(draft!.forkedAtHash).toBe("sha256:server");
    expect(draft!.id).toBeUndefined();
    expect(draft!.baseHash).toBeUndefined();
    expect(draft!.title).toBe("Fork of Server title");
    // Members come from server state (2), NOT local draft (5).
    expect(draft!.members).toHaveLength(2);
    for (const m of draft!.members) expect(m.id).toBeUndefined();

    expect(useAppState.getState().pendingConflict).toBeNull();
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/new");
  });

  it("Esc closes the modal without acting (draft preserved)", async () => {
    const user = userEvent.setup();
    const server = buildComparison({ id: 42 });
    seedDraft({ title: "Local title", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    expect(screen.getByTestId("conflict-modal")).toBeInTheDocument();
    await user.keyboard("{Escape}");
    expect(useAppState.getState().pendingConflict).toBeNull();
    // Draft preserved — Esc is non-destructive.
    expect(useAppState.getState().activeDraft).not.toBeNull();
    expect(useAppState.getState().activeDraft!.title).toBe("Local title");
    // Path unchanged.
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42/edit");
  });

  it("clicking the backdrop closes without acting (draft preserved)", async () => {
    const user = userEvent.setup();
    const server = buildComparison({ id: 42 });
    seedDraft({ title: "Local title", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    // The backdrop is the modal's outer wrapper — click it (not the inner
    // dialog) to trigger the close.
    const backdrop = screen.getByTestId("conflict-modal");
    await user.click(backdrop);
    expect(useAppState.getState().pendingConflict).toBeNull();
    expect(useAppState.getState().activeDraft).not.toBeNull();
  });

  it("focus trap cycles Tab among the three action buttons", async () => {
    const user = userEvent.setup();
    const server = buildComparison({ id: 42 });
    seedDraft({ title: "Local title", members: 3 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    const discard   = screen.getByTestId("conflict-discard");
    const fork      = screen.getByTestId("conflict-fork");
    const overwrite = screen.getByTestId("conflict-overwrite");
    discard.focus();
    expect(document.activeElement).toBe(discard);
    await user.tab();
    expect(document.activeElement).toBe(fork);
    await user.tab();
    expect(document.activeElement).toBe(overwrite);
    // Tab past last → wraps to first.
    await user.tab();
    expect(document.activeElement).toBe(discard);
  });

  it("renders without exploding when local draft is null (e.g., user navigated away)", () => {
    const server = buildComparison({ id: 42, members: 3 });
    // No draft seeded — global modal must still render the server side.
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    expect(screen.getByTestId("conflict-modal")).toBeInTheDocument();
    expect(screen.getByTestId("conflict-panel-server-members")).toHaveTextContent("3");
    expect(screen.getByTestId("conflict-panel-local-members")).toHaveTextContent("0");
  });
});
