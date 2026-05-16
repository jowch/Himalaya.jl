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
import { QueryClient, QueryClientProvider, useQueryClient } from "@tanstack/react-query";
import { useEffect, type ReactNode } from "react";
import { ConflictModal } from "../src/components/ConflictModal";
import { useAppState, LS_KEY } from "../src/state";
import { ConflictError } from "../src/api";
import type { Comparison } from "../src/api";
import { queryKeys } from "../src/queries";
import { COMPARE_DRAFT_KEY } from "../src/lib/comparison/draft";
import {
  attachConflictBridge, _resetConflictBridgeForTest,
} from "../src/lib/queue/conflictBridge";

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
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: null,
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

/**
 * Mounts the App-level MutationCache → Zustand bridge for the duration of
 * the render, mirroring what App.tsx does in production. Without this,
 * `useSaveComparison`'s 409 throws never reach `pendingConflict` (the
 * bridge is no longer in the hook). See `lib/queue/conflictBridge.ts`.
 */
function ConflictBridgeMount(): null {
  const queryClient = useQueryClient();
  const setPendingConflict = useAppState((s) => s.setPendingConflict);
  useEffect(() => {
    return attachConflictBridge(queryClient.getMutationCache(), setPendingConflict);
  }, [queryClient, setPendingConflict]);
  return null;
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
              <ConflictBridgeMount />
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

/**
 * Pre-seed all four cache keys (exposure, peaks, indices, groups) for the
 * given draft-member exposure ids so `buildOverwritePayload`'s cold-cache
 * prefetch (#74) is a no-op. Use this in tests that exercise Overwrite but
 * aren't asserting prefetch behavior — without it, the prefetch would issue
 * extra `fetch` calls that disturb the captured-request count those tests
 * pin.
 *
 * `seedDraft` numbers exposure ids starting at 300; pass that range here
 * (or use `seedWarmCachesForSeedDraft(qc, n)`).
 */
function seedWarmCaches(qc: QueryClient, exposureIds: number[]): void {
  for (const id of exposureIds) {
    qc.setQueryData(queryKeys.exposure(id), {
      id, sample_id: 1, exposure_type: "simple", filename: `exp${id}`,
    });
    qc.setQueryData(queryKeys.peaks(id),  []);
    qc.setQueryData(queryKeys.indices(id), []);
    qc.setQueryData(queryKeys.groups(id), []);
  }
}

function seedWarmCachesForSeedDraft(qc: QueryClient, members: number): void {
  seedWarmCaches(qc, Array.from({ length: members }, (_, i) => 300 + i));
}

function seedDraft(opts: {
  title?: string;
  members?: number;
  id?: number;
  baseHash?: string;
  viewGroupingMode?: "bySample" | "byPhase" | "distinct";
  viewShowPeakTicks?: boolean;
  viewShowPeakLabels?: boolean;
}) {
  const {
    title = "Local draft title", members = 3, id = 42, baseHash = "sha256:local",
    viewGroupingMode = undefined, viewShowPeakTicks = undefined, viewShowPeakLabels = undefined,
  } = opts;
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
      viewGroupingMode,
      viewShowPeakTicks,
      viewShowPeakLabels,
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
  // Reset the module-scoped bridged-mutationId set so a stale cache entry
  // from a previous test doesn't suppress legitimate bridging in the next.
  _resetConflictBridgeForTest();
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
    seedDraft({
      title: "Local title", members: 1, id: 42, baseHash: "sha256:stale",
      viewGroupingMode: "bySample", viewShowPeakTicks: true, viewShowPeakLabels: false,
    });

    const captured: { url: string; body: unknown }[] = [];
    globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = typeof input === "string" ? input : String(input);
      const body = typeof init?.body === "string" ? JSON.parse(init.body) : undefined;
      captured.push({ url, body });
      return new Response(JSON.stringify(buildComparison({
        id: 42, hash: "sha256:after-overwrite", title: "Local title",
      })), { status: 200, headers: { "Content-Type": "application/json" } });
    }) as typeof fetch;

    const { qc } = renderModal();
    // Pre-warm caches so #74's prefetch is a no-op — this test pins the
    // captured-request count and a prefetch would inflate it.
    seedWarmCachesForSeedDraft(qc, 1);
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
      view_grouping_mode:    "bySample",
      view_show_peak_ticks:  true,
      view_show_peak_labels: false,
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

    const { qc } = renderModal();
    // Pre-warm caches so #74's prefetch is a no-op (it would otherwise hit
    // the same fake fetch and trigger a 409 cascade before the save fires).
    seedWarmCachesForSeedDraft(qc, 1);
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
    // Members come from server state (2), NOT local draft (5).
    expect(draft!.members).toHaveLength(2);
    for (const m of draft!.members) expect(m.id).toBeUndefined();

    expect(useAppState.getState().pendingConflict).toBeNull();
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/new");
  });

  // Issue #58: Fork-to-new must preserve the user's drafted title rather
  // than silently overwriting it with `Fork of <server-title>`. The fork is
  // offered as the SAFE conflict resolution; discarding the user's title
  // edit undermines that promise.
  it("Fork preserves the user's drafted title (issue #58)", async () => {
    const user = userEvent.setup();
    const server = buildComparison({
      id: 42, hash: "sha256:server", title: "Reference traces (A3)", members: 2,
    });
    seedDraft({ title: "Reference traces (B3)", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    await user.click(screen.getByTestId("conflict-fork"));

    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    expect(draft!.title).toBe("Reference traces (B3)");
  });

  it("Fork falls back to `Fork of <parent>` when the user's drafted title is empty", async () => {
    const user = userEvent.setup();
    const server = buildComparison({
      id: 42, hash: "sha256:server", title: "Server title", members: 2,
    });
    // Whitespace-only counts as empty — the helper trims before checking.
    seedDraft({ title: "   ", members: 5 });
    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", server),
      );
    });
    await user.click(screen.getByTestId("conflict-fork"));

    const draft = useAppState.getState().activeDraft;
    expect(draft!.title).toBe("Fork of Server title");
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

  it("Overwrite is debounced against double-click: two synchronous clicks fire ONE mutate", async () => {
    // queue-reviewer Fix 3: the `disabled={save.isPending}` flag flips
    // async, so two synchronous clicks both pass the disabled check and
    // both call `save.mutate(...)`. Each mints a fresh client_op_id so
    // the queue's idempotency layer doesn't dedupe them — they race
    // against the same `expected_content_hash` and the second hits a 409
    // with the user's own just-committed state showing as "server".
    //
    // Guard via `overwriteInFlightRef` at the top of `handleOverwrite`:
    // the second click reads the ref synchronously and bails before
    // mutate runs.
    const server = buildComparison({ id: 42, hash: "sha256:server-v1" });
    seedDraft({ title: "Local title", members: 1, id: 42, baseHash: "sha256:stale" });

    const captured: { url: string; body: unknown }[] = [];
    let resolveResponse: ((r: Response) => void) | undefined;
    globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = typeof input === "string" ? input : String(input);
      const body = typeof init?.body === "string" ? JSON.parse(init.body) : undefined;
      captured.push({ url, body });
      // Hold the response open so a second click happens while the first
      // mutation is still pending — closing the window the bug would have
      // exploited (the ref must hold past mutate() return, not just
      // until the response lands).
      return new Promise<Response>((resolve) => { resolveResponse = resolve; });
    }) as typeof fetch;

    const { qc } = renderModal();
    // Pre-warm caches so #74's prefetch is a no-op — this test pins the
    // captured-request count and the prefetch would inflate it.
    seedWarmCachesForSeedDraft(qc, 1);
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server-v1", server),
      );
    });
    const overwrite = screen.getByTestId("conflict-overwrite");

    // Two synchronous clicks. Don't use userEvent (which awaits) — we
    // want the back-to-back-fire path the user can hit with a fast
    // double-click before React re-renders the disabled state.
    act(() => {
      overwrite.click();
      overwrite.click();
    });

    // Settle React effects so the mutation lifecycle has had a chance to
    // run.
    await waitFor(() => expect(captured.length).toBeGreaterThanOrEqual(1));
    // Pin: only ONE request hit the network. Without the guard there
    // would be two captured submits with two distinct client_op_ids.
    expect(captured).toHaveLength(1);
    expect(captured[0]!.url).toBe("/api/comparisons/42/submit");

    // Cleanup so the held promise doesn't leak past the test.
    if (resolveResponse) {
      resolveResponse(new Response(JSON.stringify(buildComparison({
        id: 42, hash: "sha256:after-overwrite",
      })), { status: 200, headers: { "Content-Type": "application/json" } }));
    }
  });

  it("Overwrite prefetches cold member caches before computing snapshots (#74)", async () => {
    // `buildOverwritePayload` calls `computeMemberSnapshot`, which reads four
    // cache keys (exposure, peaks, indices, groups) per member. If any are
    // cold, the snapshot lands `analysis_inputs_hash = ""` and the server
    // marks the member stale on first view fold — silently regressing the
    // #49 prefetch fix that `handleSave` carries. Today the bug is masked
    // because Overwrite is reached via Save→409 (caches already warm); this
    // test exercises the cold-cache scenario to lock the fix.
    const user = userEvent.setup();
    const server = buildComparison({ id: 42, hash: "sha256:server-v1" });
    // Seed a draft member with exposure_id = 999 — NOT pre-seeded in the
    // QueryClient cache. All four cache keys for 999 must be cold so the
    // prefetch actually fires.
    useAppState.setState({
      activeDraft: {
        id: 42,
        baseHash: "sha256:stale",
        title: "Local title",
        description: "",
        members: [{
          id: undefined,
          exposure_id: 999,
          display_order: 0,
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
        }],
        forkedFromId: undefined,
        forkedAtHash: undefined,
        viewGroupingMode: undefined,
        viewShowPeakTicks: undefined,
        viewShowPeakLabels: undefined,
      },
    });

    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation((async (
      input: RequestInfo | URL,
    ) => {
      const url = typeof input === "string" ? input : String(input);
      // Prefetch endpoints — return minimal shapes that satisfy the typed
      // fetchers in api.ts (lists return arrays; exposure returns a row).
      if (/\/api\/exposures\/999(?:\?|$)/.test(url)) {
        return new Response(JSON.stringify({
          id: 999, sample_id: 1, exposure_type: "simple", filename: "x",
          analysis_inputs_hash: "sha256:exposure999",
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url.includes("/api/exposures/999/peaks")
          || url.includes("/api/exposures/999/indices")
          || url.includes("/api/exposures/999/groups")) {
        return new Response("[]", {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      // Save mutation
      return new Response(JSON.stringify(buildComparison({
        id: 42, hash: "sha256:after-overwrite", title: "Local title",
      })), { status: 200, headers: { "Content-Type": "application/json" } });
    }) as typeof fetch);

    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server-v1", server),
      );
    });

    await user.click(screen.getByTestId("conflict-overwrite"));

    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : String(c[0]),
      );
      expect(urls.some((u) => /\/api\/exposures\/999(?:\?|$)/.test(u))).toBe(true);
      expect(urls.some((u) => u.includes("/api/exposures/999/peaks"))).toBe(true);
      expect(urls.some((u) => u.includes("/api/exposures/999/indices"))).toBe(true);
      expect(urls.some((u) => u.includes("/api/exposures/999/groups"))).toBe(true);
    });

    // Sanity: the save mutation must run AFTER the prefetch — the Overwrite
    // POST should appear in the fetch log too.
    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : String(c[0]),
      );
      expect(urls.some((u) => u.includes("/api/comparisons/42/submit"))).toBe(true);
    });
  });

  it("Overwrite surfaces a toast when prefetch fails (#92 review)", async () => {
    // PR #92 review point — the catch block previously swallowed the
    // prefetch error silently. Verify the user gets feedback now.
    const user = userEvent.setup();
    const server = buildComparison({ id: 42, hash: "sha256:server-v1" });
    useAppState.setState({
      activeDraft: {
        id: 42, baseHash: "sha256:stale", title: "Local title", description: "",
        members: [{
          id: undefined, exposure_id: 999, display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: undefined, label_override: undefined,
          q_window_min: undefined, q_window_max: undefined, peak_display: undefined,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "" },
        }],
        forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined,
      },
    });

    // Make the prefetch fail. Any of the four endpoints throwing causes
    // Promise.all to reject, which our catch block handles.
    vi.spyOn(global, "fetch").mockImplementation((async () => {
      throw new Error("network down");
    }) as typeof fetch);

    // Capture the toast surface — toast.ts publishes via setToastImpl.
    const { setToastImpl } = await import("../src/lib/toast");
    const toastSpy = vi.fn();
    setToastImpl(toastSpy);

    renderModal();
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server-v1", server),
      );
    });

    await user.click(screen.getByTestId("conflict-overwrite"));

    await waitFor(() => {
      expect(toastSpy).toHaveBeenCalledWith(
        expect.stringMatching(/refresh comparison data/i),
        "error",
      );
    });

    // Restore default impl for downstream tests.
    setToastImpl(null);
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
