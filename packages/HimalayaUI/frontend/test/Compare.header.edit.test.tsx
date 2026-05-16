/**
 * Compare header (edit mode) — Compare UX C-13.
 *
 * Verifies that ComparePageEdit's edit-mode header is now built from the new
 * leaf components (`CompareTitleStrip`, `CompareStatusSurface`,
 * `CompareToolbar`) and that the `SavePill` is mounted in the toolbar's
 * `saveControl` slot — replacing the legacy `compare-edit-title` input and
 * the `comparison-save` / `comparison-cancel` / `comparison-discard` triplet.
 *
 * Mirrors the fixture scaffolding of `test/Compare.header.review.test.tsx`
 * and `test/ComparePageEdit.test.tsx`.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Compare } from "../src/pages/Compare";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type { Peak, Exposure } from "../src/api";

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

function renderEdit(opts: { qc: QueryClient; initialPath: string }) {
  return render(
    <QueryClientProvider client={opts.qc}>
      <MemoryRouter initialEntries={[opts.initialPath]}>
        <Routes>
          <Route path="/experiments/:eid/compare" element={<div data-testid="list-page" />} />
          <Route path="/experiments/:eid/compare/new" element={<Compare />} />
          <Route path="/experiments/:eid/compare/:id" element={<div data-testid="review-page" />} />
          <Route path="/experiments/:eid/compare/:id/edit" element={<Compare />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare edit header — Compare UX C-13", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    sessionStorage.clear();
    localStorage.clear();
    useAppState.setState({ activeDraft: null, username: "alice" });
  });

  it("renders the new title strip + toolbar (not the legacy title input)", async () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("My comparison");
    useAppState.getState().addMember(100, qc);

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });

    expect(await screen.findByTestId("compare-title-strip")).toBeInTheDocument();
    expect(await screen.findByTestId("compare-toolbar")).toBeInTheDocument();

    // Legacy title input is gone from the edit header.
    expect(screen.queryByTestId("compare-edit-title")).toBeNull();
    expect(screen.queryByTestId("compare-edit-description")).toBeNull();
  });

  it("renders Save pill in place of Save/Cancel/Discard triplet", async () => {
    const qc = makeQc();
    seedExposure(qc, 100);
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("My comparison");
    useAppState.getState().addMember(100, qc);

    renderEdit({ qc, initialPath: "/experiments/7/compare/new" });

    expect(await screen.findByTestId("save-pill")).toBeInTheDocument();
    expect(screen.queryByTestId("comparison-save")).toBeNull();
    expect(screen.queryByTestId("comparison-cancel")).toBeNull();
    expect(screen.queryByTestId("comparison-discard")).toBeNull();
  });
});
