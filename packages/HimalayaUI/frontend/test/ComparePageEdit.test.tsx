/**
 * ComparePageEdit (Plan §Phase 4, Task 4.4).
 *
 * Asserts:
 * - Save button disabled when draft has 0 members.
 * - Save calls saveComparison with the right payload (id absent ⇒ create,
 *   id present ⇒ submit with expected_content_hash from baseHash).
 * - Cancel navigation: from /new ⇒ /experiments/:eid/compare; from
 *   /:id/edit ⇒ /experiments/:eid/compare/:id.
 * - Discard clears draft and navigates back to list.
 *
 * The save flow uses the real `useSaveComparison` hook. We mock the
 * underlying `saveComparison` via `global.fetch` per the JSDOM interceptor
 * pattern in CLAUDE.md.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparePageEdit } from "../src/pages/ComparePageEdit";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type {
  Peak, Exposure, Comparison,
} from "../src/api";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function PathProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="path-probe">{loc.pathname}</div>;
}

function renderEdit(opts: {
  qc: QueryClient;
  initialPath: string;
}) {
  return render(
    <QueryClientProvider client={opts.qc}>
      <MemoryRouter initialEntries={[opts.initialPath]}>
        <Routes>
          <Route path="/experiments/:eid/compare" element={<><PathProbe /><div data-testid="list-page" /></>} />
          <Route path="/experiments/:eid/compare/new" element={<><PathProbe /><ComparePageEdit /></>} />
          <Route path="/experiments/:eid/compare/:id" element={<><PathProbe /><div data-testid="review-page" /></>} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<><PathProbe /><ComparePageEdit /></>} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function seedExposure(qc: QueryClient, exposureId: number): void {
  const peaks: Peak[] = [
    { id: 1, exposure_id: exposureId, q: 0.10, intensity: 1.0, sharpness: 0.5, prominence: null, source: "auto", excluded: false },
  ];
  const exposure: Exposure = {
    id: exposureId,
    sample_id: 1,
    filename: "x.dat",
    kind: "file",
    selected: false,
    status: "accepted",
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: "tr",
    analysis_inputs_hash: "abcd",
  };
  qc.setQueryData(queryKeys.peaks(exposureId), peaks);
  qc.setQueryData(queryKeys.indices(exposureId), []);
  qc.setQueryData(queryKeys.groups(exposureId), []);
  qc.setQueryData(queryKeys.exposure(exposureId), exposure);
  // Seed trace cache so `useMemberTracesLoading` resolves false and the
  // `compare-edit-plot` skeleton wraps lift, exposing the gutter testid.
  qc.setQueryData(["exposure", exposureId, "trace"] as const, {
    q: [0.1, 0.2], I: [1.0, 0.5], sigma: [0.01, 0.01],
  });
}

beforeEach(() => {
  vi.restoreAllMocks();
  sessionStorage.clear();
  localStorage.clear();
  useAppState.setState({ activeDraft: null, username: "alice" });
});

describe("ComparePageEdit", () => {
  it("Save button is disabled when the draft has zero members", () => {
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    expect(screen.getByTestId("comparison-save")).toBeDisabled();
  });

  it("Save with members posts to /api/comparisons (create flow) and navigates to review", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("My comparison");
    useAppState.getState().addMember(100, qc);

    const created: Comparison = {
      id: 42,
      title: "My comparison",
      description: null,
      content_hash: "h-new",
      created_by: 1,
      created_at: "2026-05-01T00:00:00Z",
      updated_at: "2026-05-01T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      view_grouping_mode: null,
      view_show_peak_ticks: null,
      view_show_peak_labels: null,
      last_event_at: null,
      members: [],
    };
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/comparisons" && init?.method === "POST") {
        return new Response(JSON.stringify(created), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response("not found", { status: 404 });
    });

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    await user.click(screen.getByTestId("comparison-save"));

    await waitFor(() => {
      expect(fetchSpy).toHaveBeenCalled();
    });
    // First call: POST /api/comparisons (create)
    const firstCall = fetchSpy.mock.calls.find(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return url === "/api/comparisons" && (init as RequestInit | undefined)?.method === "POST";
    });
    expect(firstCall).toBeDefined();
    const body = JSON.parse((firstCall![1] as RequestInit).body as string);
    expect(body.title).toBe("My comparison");
    expect(body.members).toHaveLength(1);
    expect(body.expected_content_hash).toBeUndefined(); // create path

    // Successful save → navigate to review URL
    await waitFor(() => {
      expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42");
    });
    // Draft cleared
    expect(useAppState.getState().activeDraft).toBeNull();
  });

  it("Save on existing comparison posts to /:id/submit with expected_content_hash from baseHash", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    // Inject a draft that mimics post-loadDraftFromComparison state.
    useAppState.setState({
      activeDraft: {
        id: 42,
        baseHash: "h-existing",
        title: "Existing",
        description: "",
        members: [{
          id: 7,
          exposure_id: 100,
          display_order: 0,
          band_height: 1.0,
          y_offset: 0.0,
          normalization: "none",
          color_override: undefined,
          label_override: undefined,
          q_window_min: undefined,
          q_window_max: undefined,
          peak_display: undefined,
          snapshot: undefined,
        }],
        forkedFromId: undefined,
        forkedAtHash: undefined,
        viewGroupingMode: undefined,
        viewShowPeakTicks: undefined,
        viewShowPeakLabels: undefined,
      },
    });
    const updated: Comparison = {
      id: 42, title: "Existing", description: null, content_hash: "h-new",
      created_by: 1, created_at: null, updated_at: "2026-05-02T00:00:00Z",
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
      last_event_at: null,
      members: [],
    };
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/comparisons/42/submit" && init?.method === "POST") {
        return new Response(JSON.stringify(updated), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response("not found", { status: 404 });
    });

    renderEdit({ qc, initialPath: "/experiments/7/compare/42/edit" });
    await user.click(screen.getByTestId("comparison-save"));

    await waitFor(() => {
      expect(fetchSpy).toHaveBeenCalled();
    });
    const submitCall = fetchSpy.mock.calls.find(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return url === "/api/comparisons/42/submit" && (init as RequestInit | undefined)?.method === "POST";
    });
    expect(submitCall).toBeDefined();
    const body = JSON.parse((submitCall![1] as RequestInit).body as string);
    expect(body.title).toBe("Existing");
    expect(body.expected_content_hash).toBe("h-existing");
    expect(body.members).toHaveLength(1);

    await waitFor(() => {
      expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42");
    });
  });

  it("Cancel from /new returns to /experiments/:eid/compare", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    await user.click(screen.getByTestId("comparison-cancel"));
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare");
  });

  it("Cancel from /:id/edit returns to /experiments/:eid/compare/:id", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    useAppState.setState({
      activeDraft: {
        id: 42,
        baseHash: "h",
        title: "x",
        description: "",
        members: [],
        forkedFromId: undefined,
        forkedAtHash: undefined,
        viewGroupingMode: undefined,
        viewShowPeakTicks: undefined,
        viewShowPeakLabels: undefined,
      },
    });
    renderEdit({ qc, initialPath: "/experiments/7/compare/42/edit" });
    await user.click(screen.getByTestId("comparison-cancel"));
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare/42");
  });

  it("Discard clears draft and navigates to list", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("Sticky");
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    await user.click(screen.getByTestId("comparison-discard"));
    expect(useAppState.getState().activeDraft).toBeNull();
    expect(screen.getByTestId("path-probe")).toHaveTextContent("/experiments/7/compare");
  });

  it("right slot hosts ComparisonPickerPanel in edit mode", async () => {
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    expect(await screen.findByTestId("comparison-picker-panel")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-edit-right-hint")).toBeNull();
  });

  it("empty-state '+ Add traces' button focuses the panel's search input", async () => {
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    // Wait for panel to mount.
    await screen.findByTestId("comparison-picker-panel");
    fireEvent.click(screen.getByTestId("compare-edit-add-traces"));
    expect(screen.getByTestId("comparison-picker-search")).toHaveFocus();
  });

  // ── Phase 7 wiring ───────────────────────────────────────────────────────

  it("Reset heights button is rendered in the edit-mode header", () => {
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    expect(screen.getByTestId("compare-edit-reset-heights")).toBeInTheDocument();
  });

  it("Reset heights button restores all band_heights to 1.0", () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(100, qc);
    // Inflate band_height to verify the reset.
    useAppState.getState().updateMember(0, { band_height: 2.5 });
    expect(useAppState.getState().activeDraft!.members[0]!.band_height).toBe(2.5);

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    fireEvent.click(screen.getByTestId("compare-edit-reset-heights"));
    expect(useAppState.getState().activeDraft!.members[0]!.band_height).toBe(1);
  });

  it("metadata gutter mounts with the draft members in edit mode", () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(100, qc);
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    expect(screen.getByTestId("compare-edit-gutter")).toBeInTheDocument();
    // Member meta row mounted via the gutter
    expect(screen.getByTestId("member-meta-row")).toBeInTheDocument();
  });

  // ─── Keyboard shortcuts (Phase 13, Task 13.4) ─────────────────────────────

  it("Escape cancels the edit (returns to /:id review for edit-existing)", () => {
    const qc = makeQc();
    // Seed a draft tied to id=42 to exercise the cancel-from-edit path.
    useAppState.setState({
      activeDraft: {
        id: 42, baseHash: "h", title: "T", description: "",
        members: [],
        forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined,
      },
    });
    renderEdit({ qc, initialPath: "/experiments/7/compare/42/edit" });
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.getByTestId("path-probe")).toHaveTextContent(
      "/experiments/7/compare/42",
    );
  });

  it("Cmd+Enter triggers save when draft has members", async () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("kbd save");
    useAppState.getState().addMember(100, qc);

    const created: Comparison = {
      id: 99,
      title: "kbd save",
      description: null,
      content_hash: "h-new",
      created_by: 1,
      created_at: "2026-05-06T00:00:00Z",
      updated_at: "2026-05-06T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      view_grouping_mode: null,
      view_show_peak_ticks: null,
      view_show_peak_labels: null,
      last_event_at: null,
      members: [],
    };
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/comparisons" && init?.method === "POST") {
        return new Response(JSON.stringify(created), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response("not found", { status: 404 });
    });
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    fireEvent.keyDown(window, { key: "Enter", metaKey: true });
    await waitFor(() => {
      const calls = fetchSpy.mock.calls.map((c) =>
        [typeof c[0] === "string" ? c[0] : String(c[0]),
         (c[1] as RequestInit | undefined)?.method] as const,
      );
      expect(calls).toContainEqual(["/api/comparisons", "POST"]);
    });
  });

  // ── Regression: ResizeObserver re-attach when Skeleton swaps in (issue #51-edit)
  // Mirrors the review-mode #51 fix in ComparePage.test.tsx. The page-level
  // ResizeObserver gates `panelHeight`, which feeds MemberMetaGutter's
  // outer style.height. With the buggy `[plotMembers.length]` deps, the
  // effect runs once with `plotColRef.current === null` (Skeleton fallback
  // is up), bails at the early-out, and never re-fires when Skeleton lifts —
  // so the gutter rows stack at top:0 with height 0.
  it("attaches ResizeObserver after Skeleton lifts (tracesLoading flips to false)", async () => {
    // Seed peaks/exposure/etc. but NOT the trace cache — trace data is the
    // gate for `tracesLoading`. Initial mount has Skeleton up (plotColRef
    // null); when the async trace fetch lands, Skeleton lifts and the
    // observer must re-attach. Mirrors the review-mode #51 test pattern.
    const qc = makeQc();
    const peaks: Peak[] = [
      { id: 1, exposure_id: 200, q: 0.10, intensity: 1.0, sharpness: 0.5, prominence: null, source: "auto", excluded: false },
    ];
    const exposure: Exposure = {
      id: 200, sample_id: 1, filename: "x.dat", kind: "file", selected: false,
      status: "accepted", image_path: null, image_version: "", tags: [],
      sources: [], trace_hash: "tr", analysis_inputs_hash: "abcd",
    };
    qc.setQueryData(queryKeys.peaks(200), peaks);
    qc.setQueryData(queryKeys.indices(200), []);
    qc.setQueryData(queryKeys.groups(200), []);
    qc.setQueryData(queryKeys.exposure(200), exposure);

    useAppState.setState({
      activeDraft: {
        id: undefined, baseHash: undefined, title: "T", description: "",
        members: [{
          id: 1, exposure_id: 200, display_order: 0,
          band_height: 1, y_offset: 0, normalization: "max",
          color_override: undefined, label_override: undefined,
          q_window_min: undefined, q_window_max: undefined,
          peak_display: undefined,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "abcd" },
        }],
        forkedFromId: undefined, forkedAtHash: undefined,
        viewGroupingMode: undefined, viewShowPeakTicks: undefined, viewShowPeakLabels: undefined,
      },
    });

    // Async-resolve the trace fetch so Skeleton starts up, then lifts.
    vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/exposures/200/trace") {
        return new Response(
          JSON.stringify({ q: [0.1, 0.2], I: [1.0, 0.5], sigma: [0.01, 0.01] }),
          { status: 200, headers: { "Content-Type": "application/json" } },
        );
      }
      return new Response("not found", { status: 404 });
    });

    // ResizeObserver stub that fires its callback immediately on observe()
    // and surfaces the observed element's clientHeight.
    vi.stubGlobal("ResizeObserver", vi.fn((cb: ResizeObserverCallback) => {
      const inst = {
        observe: vi.fn((el: Element) => {
          Object.defineProperty(el, "clientHeight", { value: 400, configurable: true });
          cb(
            [{ target: el, contentRect: { height: 400 } } as unknown as ResizeObserverEntry],
            inst as unknown as ResizeObserver,
          );
        }),
        disconnect: vi.fn(),
      };
      return inst;
    }));

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });

    // After Skeleton lifts, the gutter container exists and the page-level
    // observer must have fired with a non-zero height. The gutter's inline
    // style.height is driven by `panelHeight` — pinning to it asserts that
    // the page-level observer attached, NOT MultiTracePlot's own observer.
    await waitFor(() => {
      const gutter = screen.getByTestId("member-meta-gutter");
      expect(gutter.style.height).not.toBe("0px");
    });
  });

  it("Cmd+Enter is a no-op when draft has zero members", () => {
    const qc = makeQc();
    useAppState.getState().startNewDraft();
    const fetchSpy = vi.spyOn(global, "fetch")
      .mockResolvedValue(new Response("[]", { status: 200,
        headers: { "Content-Type": "application/json" } }));
    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    fireEvent.keyDown(window, { key: "Enter", metaKey: true });
    // No POST /api/comparisons call. (Other reads — pins, etc — may fire,
    // but Save itself must not.)
    const saveCalls = fetchSpy.mock.calls.filter(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return url === "/api/comparisons" && (init as RequestInit | undefined)?.method === "POST";
    });
    expect(saveCalls).toHaveLength(0);
  });

  // ── PR2: right slot is always the inline ComparisonPickerPanel ───────────

  // ── Regression: cold-exposure snapshot prefetch (issue #49) ────────────────

  it("Save prefetches exposure, peaks, indices, and groups for cold members", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    // Do NOT seed exposure 200 — it is "cold" in the cache.
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("Cold member test");
    useAppState.getState().addMember(200, qc);

    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
      const url = typeof input === "string" ? input : String(input);
      if (url === "/api/exposures/200" && init?.method === "GET") {
        return new Response(JSON.stringify({
          id: 200, sample_id: 1, filename: "cold.dat", kind: "file",
          selected: false, status: "accepted", image_path: null,
          image_version: "", tags: [], sources: [], trace_hash: "tr",
          analysis_inputs_hash: "cold-hash",
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/200/peaks" && init?.method === "GET") {
        return new Response(JSON.stringify([
          { id: 1, exposure_id: 200, q: 0.10, intensity: 1.0, sharpness: 0.5, source: "auto", excluded: false },
        ]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/200/indices" && init?.method === "GET") {
        return new Response(JSON.stringify([]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/exposures/200/groups" && init?.method === "GET") {
        return new Response(JSON.stringify([]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url === "/api/comparisons" && init?.method === "POST") {
        return new Response(JSON.stringify({
          id: 77, title: "Cold member test", description: null,
          content_hash: "h-new", created_by: 1,
          created_at: "2026-05-01T00:00:00Z", updated_at: "2026-05-01T00:00:00Z",
          forked_from_id: null, forked_at_hash: null, forked_from_title: null,
          members: [],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });
    await user.click(screen.getByTestId("comparison-save"));

    await waitFor(() => {
      const calls = fetchSpy.mock.calls.map((c) =>
        [typeof c[0] === "string" ? c[0] : String(c[0]),
         (c[1] as RequestInit | undefined)?.method] as const,
      );
      expect(calls).toContainEqual(["/api/comparisons", "POST"]);
    });

    // Assert that all four cold-read GETs happened before the save POST.
    const callOrder = fetchSpy.mock.calls
      .filter(([u]) => {
        const url = typeof u === "string" ? u : String(u);
        return url.startsWith("/api/exposures/200") || url === "/api/comparisons";
      })
      .map(([u, init]) => ({
        url: typeof u === "string" ? u : String(u),
        method: (init as RequestInit | undefined)?.method ?? "GET",
      }));

    const saveIndex = callOrder.findIndex((c) => c.url === "/api/comparisons" && c.method === "POST");
    expect(saveIndex).toBeGreaterThanOrEqual(4); // exposure + peaks + indices + groups

    // Assert the submitted snapshot carries the warmed data.
    const saveCall = fetchSpy.mock.calls.find(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return url === "/api/comparisons" && (init as RequestInit | undefined)?.method === "POST";
    });
    const body = JSON.parse((saveCall![1] as RequestInit).body as string);
    expect(body.members).toHaveLength(1);
    expect(body.members[0].snapshot.analysis_inputs_hash).toBe("cold-hash");
    expect(body.members[0].snapshot.effective_peaks).toHaveLength(1);
  });
});

// ── Issue #69: edit-mode cold-load cache hydration parity ─────────────────────
//
// Mirror of the #61/#52 review-mode regression in test/ComparePage.test.tsx.
// Edit mode is a separate page component and used to skip the per-member
// exposure + sample subscriptions, so deep-linking
// `/experiments/:eid/compare/:id/edit` cold (no Inspect visit warming the
// cache, no picker pre-warm) produced gray (ORPHAN_FALLBACK) traces and
// gutter rows labeled `Exposure #N`. Fix: subscribe to useMemberExposures
// and useMemberSamples in ComparePageEdit and pipe the resolved labels into
// MemberMetaGutter via the shared resolver in lib/comparison/labels.ts.

describe("ComparePageEdit — cold-cache exposure + sample hydration (#69)", () => {
  beforeEach(() => {
    vi.stubGlobal("ResizeObserver", vi.fn(() => ({
      observe: vi.fn((el: Element) => {
        Object.defineProperty(el, "clientHeight", { value: 400, configurable: true });
      }),
      disconnect: vi.fn(),
    })));
    // Reset Zustand draft so a previous test doesn't leak an in-progress
    // draft past the URL hydration guard in ComparePageEdit (the
    // `loadDraftFromComparison` call no-ops if the active draft already
    // matches the comparison id).
    useAppState.getState().discardDraft();
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
      // Peaks/indices/groups are read by computeMemberSnapshot at submit
      // time, not on mount — return empty 200s so the load path stays clean.
      if (/^\/api\/exposures\/\d+\/(peaks|indices|groups)$/.test(url)) {
        return new Response(JSON.stringify([]), {
          status: 200, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response("not found", { status: 404 });
    }) as unknown as typeof fetch;
  }

  it("hydrates the per-member exposure cache on edit-mode cold-load (#69)", async () => {
    vi.spyOn(global, "fetch").mockImplementation(makeFetchMock());
    const qc = makeQc();

    renderEdit({ qc, initialPath: "/experiments/7/compare/42/edit" });

    // Wait for the comparison fetch to resolve and the draft to hydrate.
    await waitFor(() => {
      expect(screen.getByTestId("compare-page-edit")).toHaveAttribute(
        "data-comparison-id", "42",
      );
    });

    // Bug: ComparePageEdit never subscribed to per-member exposure rows in
    // cold-load (the picker pre-warms via #49, but deep-linking
    // /:id/edit bypasses the picker). Without an explicit subscription the
    // exposure cache stays undefined indefinitely.
    await waitFor(() => {
      expect(qc.getQueryData(queryKeys.exposure(100))).toBeDefined();
      expect(qc.getQueryData(queryKeys.exposure(200))).toBeDefined();
    });
  });

  it("renders gutter labels with sample label + filename, not 'Exposure #N' (#69)", async () => {
    vi.spyOn(global, "fetch").mockImplementation(makeFetchMock());
    const qc = makeQc();

    renderEdit({ qc, initialPath: "/experiments/7/compare/42/edit" });

    await waitFor(() => {
      expect(screen.getByTestId("compare-page-edit")).toHaveAttribute(
        "data-comparison-id", "42",
      );
    });

    // After exposure + sample data lands, the edit-mode gutter resolves
    // the same human-readable labels as review mode. The `Exposure #N`
    // form is the regressed fallback — assert it's gone.
    await waitFor(() => {
      const labels = screen.getAllByTestId("member-meta-label").map((n) => n.textContent ?? "");
      expect(labels.length).toBeGreaterThan(0);
      expect(labels.some((t) => t.includes("Lipid-A") && t.includes("run-A"))).toBe(true);
      expect(labels.some((t) => t.includes("Lipid-B") && t.includes("run-B"))).toBe(true);
      expect(labels.every((t) => !t.startsWith("Exposure #"))).toBe(true);
    });
  });
});
