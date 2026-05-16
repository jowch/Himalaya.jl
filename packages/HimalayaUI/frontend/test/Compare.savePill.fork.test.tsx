/**
 * Compare Save-as-fork flow — Compare UX C-14.
 *
 * When a NON-author opens a draft on someone else's comparison, `useCompareMode`
 * reports `editing-as-fork-of` and the SavePill copy reads "Save as fork…".
 * Clicking it must:
 *   1. prompt for a fork title (default `Copy of <parent title>`),
 *   2. morph the draft into a fork (clear id + baseHash, set
 *      forkedFromId + forkedAtHash, set the new title) via `setDraftForkOf`,
 *   3. submit via the CREATE path — POST /api/comparisons, no
 *      `expected_content_hash`, with `forked_from_id` + `forked_at_hash`.
 *
 * Mirrors the fixture scaffolding of `test/Compare.header.edit.test.tsx`
 * and `test/ComparePageEdit.test.tsx`.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparePageEdit } from "../src/pages/ComparePageEdit";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type { Peak, Exposure, Comparison, User } from "../src/api";

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

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
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
  qc.setQueryData(["exposure", exposureId, "trace"] as const, {
    q: [0.1, 0.2], I: [1.0, 0.5], sigma: [0.01, 0.01],
  });
}

/**
 * Seed the comparison + users caches so `useComparison` resolves
 * synchronously and `useCurrentUserId` resolves `bob` (id 2) — a NON-author
 * of comparison #7, which is `created_by: 1` (alice). This drives
 * `useCompareMode` to `editing-as-fork-of`.
 */
function seedComparisonAndUsers(qc: QueryClient, comp: Comparison): void {
  qc.setQueryData(queryKeys.comparison(comp.id), comp);
  qc.setQueryData(queryKeys.comparisonForks(comp.id), []);
  const users: User[] = [
    { id: 1, username: "alice", first_name: null, last_name: null },
    { id: 2, username: "bob",   first_name: null, last_name: null },
  ];
  qc.setQueryData(["users"] as const, users);
}

function renderEdit(opts: { qc: QueryClient; initialPath: string }) {
  return render(
    <QueryClientProvider client={opts.qc}>
      <MemoryRouter initialEntries={[opts.initialPath]}>
        <Routes>
          <Route path="/experiments/:eid/compare" element={<div data-testid="list-page" />} />
          <Route path="/experiments/:eid/compare/new" element={<ComparePageEdit />} />
          <Route path="/experiments/:eid/compare/:id" element={<div data-testid="review-page" />} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<ComparePageEdit />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

const PARENT: Comparison = {
  id: 7,
  title: "Alice's comparison",
  description: null,
  content_hash: "h-parent",
  created_by: 1, // alice — NOT the current user (bob)
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

/** A draft mimicking post-loadDraftFromComparison state for comparison #7. */
function forkableDraft() {
  return {
    id: 7 as number | undefined,
    baseHash: "h-parent" as string | undefined,
    title: "Alice's comparison",
    description: "",
    members: [{
      id: 7 as number | undefined,
      exposure_id: 100 as number | null,
      display_order: 0,
      band_height: 1.0,
      y_offset: 0.0,
      normalization: "none" as const,
      color_override: undefined,
      label_override: undefined,
      q_window_min: undefined,
      q_window_max: undefined,
      peak_display: undefined,
      snapshot: undefined,
    }],
    forkedFromId: undefined as number | undefined,
    forkedAtHash: undefined as string | undefined,
    viewGroupingMode: undefined,
    viewShowPeakTicks: undefined,
    viewShowPeakLabels: undefined,
  };
}

describe("Compare Save-as-fork flow — Compare UX C-14", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    sessionStorage.clear();
    localStorage.clear();
    // Current user is `bob` — a non-author of comparison #7 (created_by: 1).
    useAppState.setState({ activeDraft: null, username: "bob" });
  });

  it("submits via create path with fork lineage when user is non-author", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    seedComparisonAndUsers(qc, PARENT);
    useAppState.setState({ activeDraft: forkableDraft() });

    const promptSpy = vi
      .spyOn(window, "prompt")
      .mockReturnValue("My fork");

    const created: Comparison = {
      ...PARENT,
      id: 99,
      title: "My fork",
      content_hash: "h-fork",
      forked_from_id: 7,
      forked_at_hash: "h-parent",
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

    renderEdit({ qc, initialPath: "/experiments/7/compare/7/edit" });

    // SavePill copy for a non-author reads "Save as fork…".
    const pill = await screen.findByTestId("save-pill");
    expect(pill).toHaveTextContent(/fork/i);

    await user.click(pill);

    // The fork-title prompt fired.
    expect(promptSpy).toHaveBeenCalled();

    // The create path fired: POST /api/comparisons.
    await waitFor(() => {
      const created = fetchSpy.mock.calls.find(([u, init]) => {
        const url = typeof u === "string" ? u : String(u);
        return url === "/api/comparisons" && (init as RequestInit | undefined)?.method === "POST";
      });
      expect(created).toBeDefined();
    });
    const createCall = fetchSpy.mock.calls.find(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return url === "/api/comparisons" && (init as RequestInit | undefined)?.method === "POST";
    });
    const body = JSON.parse((createCall![1] as RequestInit).body as string);

    // Fork lineage rides the body.
    expect(body.forked_from_id).toBe(7);
    expect(body.forked_at_hash).toBe("h-parent");
    expect(body.title).toBe("My fork");
    // Create path — no id, no expected_content_hash.
    expect(body.id).toBeUndefined();
    expect(body.expected_content_hash).toBeUndefined();
    expect(body.members).toHaveLength(1);

    // No /:id/submit call was made.
    const submitCall = fetchSpy.mock.calls.find(([u]) => {
      const url = typeof u === "string" ? u : String(u);
      return /\/submit$/.test(url);
    });
    expect(submitCall).toBeUndefined();
  });

  it("does not mutate when the fork-title prompt is cancelled", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    seedExposure(qc, 100);
    seedComparisonAndUsers(qc, PARENT);
    useAppState.setState({ activeDraft: forkableDraft() });

    const promptSpy = vi.spyOn(window, "prompt").mockReturnValue(null);
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("not found", { status: 404 }),
    );

    renderEdit({ qc, initialPath: "/experiments/7/compare/7/edit" });

    const pill = await screen.findByTestId("save-pill");
    await user.click(pill);

    expect(promptSpy).toHaveBeenCalled();

    // Give any (erroneous) async submit a chance to fire, then assert none did.
    await new Promise((r) => setTimeout(r, 50));
    const saveCall = fetchSpy.mock.calls.find(([u, init]) => {
      const url = typeof u === "string" ? u : String(u);
      return (url === "/api/comparisons" || /\/submit$/.test(url))
        && (init as RequestInit | undefined)?.method === "POST";
    });
    expect(saveCall).toBeUndefined();

    // Draft untouched — still points at the parent (no morph).
    const d = useAppState.getState().activeDraft;
    expect(d?.id).toBe(7);
    expect(d?.forkedFromId).toBeUndefined();
  });
});
