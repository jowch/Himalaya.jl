/**
 * Compare empty state — Compare UX C-15 Step 4.
 *
 * When a draft has zero members the center pane renders the §7.1 empty block:
 * a "No traces yet." headline, a "+ Add traces" button, and the
 * "Or drag exposures from the picker." copy. The region is also a drop target
 * for external exposure drags.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Compare } from "../src/pages/Compare";
import { useAppState } from "../src/state";

// MultiTracePlot pulls Observable Plot — stub it to avoid jsdom weirdness.
vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => document.createElement("div")),
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

function seedEmptyDraft(): void {
  useAppState.setState({
    activeDraft: {
      id: undefined,
      baseHash: undefined,
      title: "",
      description: "",
      members: [],
      forkedFromId: undefined,
      forkedAtHash: undefined,
      viewGroupingMode: undefined,
      viewShowPeakTicks: undefined,
      viewShowPeakLabels: undefined,
    },
  });
}

function renderCompare(qc: QueryClient) {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/compare/new"]}>
        <Routes>
          <Route path="/experiments/:eid/compare/new" element={<Compare />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare empty state — Compare UX C-15", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    sessionStorage.clear();
    localStorage.clear();
    useAppState.setState({ activeDraft: null, username: "alice" });
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
  });

  it("renders the empty block when a draft has zero members", () => {
    const qc = makeQc();
    seedEmptyDraft();
    renderCompare(qc);
    expect(screen.getByTestId("compare-empty-state")).toBeInTheDocument();
    expect(screen.getByText(/No traces yet/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /Add traces/i })).toBeInTheDocument();
    expect(screen.getByText(/Or drag exposures/i)).toBeInTheDocument();
  });

  it("clicking + Add traces opens the interim picker chip", async () => {
    const qc = makeQc();
    seedEmptyDraft();
    renderCompare(qc);
    // The interim picker chip is the ComparisonPickerPanel; clicking
    // "+ Add traces" focuses its search input.
    const panel = await screen.findByTestId("comparison-picker-panel");
    expect(panel).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /Add traces/i }));
    expect(screen.getByTestId("comparison-picker-search")).toHaveFocus();
  });

  it("renders as a drop target for external exposure drags", () => {
    const qc = makeQc();
    seedEmptyDraft();
    renderCompare(qc);
    const region = screen.getByTestId("compare-empty-state");
    // A dragover whose default is prevented marks the region a valid drop
    // target. `fireEvent.dragOver` returns false when a handler called
    // preventDefault — assert the handler is wired.
    const notCancelled = fireEvent.dragOver(region, {
      dataTransfer: { dropEffect: "" },
    });
    expect(notCancelled).toBe(false);
    // The drop handler routes to the interim picker (focuses its search).
    fireEvent.drop(region, { dataTransfer: { dropEffect: "" } });
    expect(screen.getByTestId("comparison-picker-search")).toHaveFocus();
  });
});
