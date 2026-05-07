/**
 * Warm-add accelerator (Plan §Phase 5, Task 5.3).
 *
 * The Inspect page surfaces an "Add to comparison" affordance per active
 * exposure. Dropdown options for v1:
 *   - "Recent draft: <title>" — appears when activeDraft is non-null in
 *     sessionStorage. Adds the active exposure to the current draft.
 *   - "+ New comparison" — clears any existing draft, starts a new one
 *     pre-populated with the active exposure, navigates to
 *     /experiments/:eid/compare/new.
 *
 * Deferred for a follow-up (per spec): "Pick a comparison..." which would
 * open a comparison-scoped picker. The plan flags this as v1-optional.
 *
 * Cross-tab boundary test: the "Add to current draft" path mutates only
 * the same tab's Zustand draft slot (sessionStorage-scoped). A second
 * tab's draft is unchanged until submit + SSE — modeled here by direct
 * sessionStorage inspection (sessionStorage is per-tab in real browsers
 * but global in JSDOM; we assert the action only writes to sessionStorage,
 * not localStorage or any cross-tab medium).
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { WarmAddMenu } from "../src/components/WarmAddMenu";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import { COMPARE_DRAFT_KEY } from "../src/lib/comparison/draft";
import type { Exposure, Peak } from "../src/api";

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

function renderHarness(
  ui: JSX.Element,
  opts: { qc: QueryClient; initialPath: string } = {
    qc: makeQc(),
    initialPath: "/experiments/7/inspect",
  },
) {
  return render(
    <QueryClientProvider client={opts.qc}>
      <MemoryRouter initialEntries={[opts.initialPath]}>
        <PathProbe />
        <Routes>
          <Route path="/experiments/:eid/inspect" element={ui} />
          <Route path="/experiments/:eid/compare/new" element={<div data-testid="edit-page" />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function seedExposure(qc: QueryClient, exposureId: number): void {
  const peaks: Peak[] = [
    { id: 1, exposure_id: exposureId, q: 0.10, intensity: 1.0, sharpness: 0.5, source: "auto", excluded: false },
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
}

beforeEach(() => {
  sessionStorage.clear();
  localStorage.clear();
  vi.restoreAllMocks();
  useAppState.setState({
    activeDraft: null,
    username: "alice",
    activeExperimentId: 7,
  });
});

describe("WarmAddMenu", () => {
  it("does not render menu items until trigger is clicked", () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    renderHarness(
      <WarmAddMenu exposureId={100} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    expect(screen.queryByTestId("warm-add-new")).toBeNull();
  });

  it("opens menu on trigger click; '+ New comparison' navigates to /new with the exposure pre-populated", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    renderHarness(
      <WarmAddMenu exposureId={100} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    await user.click(screen.getByTestId("warm-add-trigger"));
    const newBtn = await screen.findByTestId("warm-add-new");
    await user.click(newBtn);

    // Navigated to the edit page.
    await waitFor(() => {
      expect(screen.getByTestId("path-probe"))
        .toHaveTextContent("/experiments/7/compare/new");
    });
    // Draft pre-populated.
    const d = useAppState.getState().activeDraft;
    expect(d).not.toBeNull();
    expect(d!.members).toHaveLength(1);
    expect(d!.members[0]!.exposure_id).toBe(100);
  });

  it("'Recent draft' menu item is hidden when no draft exists", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    renderHarness(
      <WarmAddMenu exposureId={100} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    await user.click(screen.getByTestId("warm-add-trigger"));
    expect(screen.queryByTestId("warm-add-recent")).toBeNull();
  });

  it("'Recent draft' shows when a draft exists; clicking adds the exposure to that draft", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    seedExposure(qc, 200);
    // Pre-existing draft with one member; adding 200 should make it two.
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "My ongoing comparison",
        description: "",
        members: [{
          id: undefined,
          exposure_id: 100,
          display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: undefined, label_override: undefined,
          q_window_min: undefined, q_window_max: undefined,
          peak_display: undefined, snapshot: undefined,
        }],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });

    renderHarness(
      <WarmAddMenu exposureId={200} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    await user.click(screen.getByTestId("warm-add-trigger"));
    const recentBtn = await screen.findByTestId("warm-add-recent");
    expect(recentBtn).toHaveTextContent(/My ongoing comparison/);
    await user.click(recentBtn);

    const d = useAppState.getState().activeDraft;
    expect(d).not.toBeNull();
    expect(d!.members).toHaveLength(2);
    expect(d!.members.map((m) => m.exposure_id)).toEqual([100, 200]);
  });

  it("'Recent draft' hides the menu item if the exposure is already in the draft", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "X",
        description: "",
        members: [{
          id: undefined,
          exposure_id: 100,
          display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: undefined, label_override: undefined,
          q_window_min: undefined, q_window_max: undefined,
          peak_display: undefined, snapshot: undefined,
        }],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });

    renderHarness(
      <WarmAddMenu exposureId={100} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    await user.click(screen.getByTestId("warm-add-trigger"));
    // Recent menu item suppressed because exposure is already a member —
    // user is given the explicit hint instead.
    const recent = screen.queryByTestId("warm-add-recent");
    expect(recent).toBeNull();
    expect(screen.getByTestId("warm-add-already-added")).toBeInTheDocument();
  });

  it("trigger is disabled when no exposure is selected", async () => {
    const qc = makeQc();
    renderHarness(
      <WarmAddMenu exposureId={undefined} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    expect((screen.getByTestId("warm-add-trigger") as HTMLButtonElement).disabled).toBe(true);
  });

  // ── Cross-tab boundary test ────────────────────────────────────────────
  it("adding to current draft writes only to sessionStorage, not localStorage", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    seedExposure(qc, 200);
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "Draft",
        description: "",
        members: [{
          id: undefined,
          exposure_id: 100,
          display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: undefined, label_override: undefined,
          q_window_min: undefined, q_window_max: undefined,
          peak_display: undefined, snapshot: undefined,
        }],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });

    // Snapshot localStorage before the action.
    const lsBefore = JSON.stringify({ ...localStorage });
    // Spy on storage setters so we can prove no cross-tab broadcast fired.
    const setItemSpy = vi.spyOn(Storage.prototype, "setItem");

    renderHarness(
      <WarmAddMenu exposureId={200} experimentId={7} />,
      { qc, initialPath: "/experiments/7/inspect" },
    );
    await user.click(screen.getByTestId("warm-add-trigger"));
    await user.click(await screen.findByTestId("warm-add-recent"));

    // Draft mutation lands in sessionStorage under COMPARE_DRAFT_KEY.
    const persisted = sessionStorage.getItem(COMPARE_DRAFT_KEY);
    expect(persisted).not.toBeNull();
    const parsed = JSON.parse(persisted!);
    expect(parsed.draft.members).toHaveLength(2);

    // Cross-tab boundary: localStorage's draft slot is untouched. We don't
    // mirror activeDraft into the persisted Zustand slice (it's not in the
    // partializer), so any localStorage changes are unrelated to the
    // warm-add action.
    const lsAfter = JSON.stringify({ ...localStorage });
    expect(lsAfter).toBe(lsBefore);

    // BroadcastChannel is intentionally NOT used in v1 (per spec) — no
    // BroadcastChannel constructor calls should have fired. JSDOM lacks
    // BroadcastChannel by default, so we just assert no localStorage
    // setItem hit a known cross-tab broadcast key.
    const broadcastWrites = setItemSpy.mock.calls.filter(([key]) =>
      typeof key === "string" && /broadcast|cross-tab/.test(key),
    );
    expect(broadcastWrites).toHaveLength(0);
  });
});
