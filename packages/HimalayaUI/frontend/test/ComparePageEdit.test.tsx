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
import { render, screen, waitFor } from "@testing-library/react";
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
      },
    });
    const updated: Comparison = {
      id: 42, title: "Existing", description: null, content_hash: "h-new",
      created_by: 1, created_at: null, updated_at: "2026-05-02T00:00:00Z",
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
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
});
