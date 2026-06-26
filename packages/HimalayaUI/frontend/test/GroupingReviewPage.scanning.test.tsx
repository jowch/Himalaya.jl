import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import type React from "react";
import { render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/pages/GroupingReviewPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useAppState } from "../src/state";

// The page reads loads (live-invalidated) + the experiment (name/ingest_status)
// + the SSE ingestInFlight. Mock the queries; drive ingestInFlight via the store.
let LOADS: Load[] = [];
let EXP = { id: 7, name: "SSRL April 2026 - 1p7m", ingest_status: "scanning" };
vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useExperiment: () => ({ data: EXP }),
  };
});

// TestShell: needed so InteractionDock renders the dock-action-confirm button.
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return <>{children}<InteractionDock /></>;
}

function wrap(node: React.ReactNode) {
  const qc = new QueryClient();
  return render(
    <QueryClientProvider client={qc}><TestShell>{node}</TestShell></QueryClientProvider>,
  );
}

const flaggedLoad: Load = {
  load_id: 6, load_index: 6, session_id: null, start_time: "10:02", end_time: "10:38",
  frame_count: 0, note: null,
  samples: [{
    sample_id: 50, name: "2-2 + LL37 1:1", slot_index: 5, grouping_source: "auto_position",
    name_source: "auto", merged_into_id: null,
    flag: { kind: "split", split_at_index: 2, jump_from: 12.4, jump_to: 48.1 }, exposures: [],
  }],
};
const cleanLoad: Load = {
  load_id: 7, load_index: 7, session_id: null, start_time: "10:40", end_time: "10:55",
  frame_count: 0, note: null,
  samples: [{
    sample_id: 60, name: "CM1 Only", slot_index: 6, grouping_source: "auto_position",
    name_source: "auto", merged_into_id: null, flag: null, exposures: [],
  }],
};

describe("GroupingReviewPage scanning surface (p1-grouping)", () => {
  beforeEach(() => {
    LOADS = [flaggedLoad, cleanLoad];
    EXP = { id: 7, name: "SSRL April 2026 - 1p7m", ingest_status: "scanning" };
    useAppState.setState({ ingestInFlight: { 7: { status: "scanning", processed: 418, total: 682 } } });
  });
  afterEach(() => useAppState.setState({ ingestInFlight: null }));

  it("scanning: shows the live header, counts, flag count, unfolding tail, disabled Confirm", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.getByTestId("grouping-scanning-header")).toBeInTheDocument();
    expect(screen.getByText(/Parsing exposures/i)).toBeInTheDocument();
    expect(screen.getByText(/418 \/ ~682/)).toBeInTheDocument();
    // One flagged sample across the two loads.
    expect(screen.getByTestId("grouping-flag-count")).toHaveTextContent("1 flag to review");
    // processed (418) < total (682) → the "unfolding…" tail shows.
    expect(screen.getByTestId("grouping-unfolding")).toBeInTheDocument();
    // Confirm is gated while scanning (dock-action-confirm rendered by InteractionDock).
    expect(screen.getByTestId("dock-action-confirm")).toBeDisabled();
    // BOTH loads are visible during a scan: the clean load is shown (collapsed,
    // "grouped cleanly") rather than filtered out the way the post-scan attn
    // filter would hide it. The flagged load's sample is expanded.
    expect(screen.getByText(/grouped cleanly/i)).toBeInTheDocument();
    expect(screen.getByText("2-2 + LL37 1:1")).toBeInTheDocument();
  });

  it("complete: no scanning header, Confirm enabled, fires onConfirm", () => {
    useAppState.setState({ ingestInFlight: null });
    EXP = { id: 7, name: "SSRL April 2026 - 1p7m", ingest_status: "complete" };
    const onConfirm = vi.fn();
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} onConfirm={onConfirm} />);
    expect(screen.queryByTestId("grouping-scanning-header")).toBeNull();
    const confirm = screen.getByTestId("dock-action-confirm");
    expect(confirm).toBeEnabled();
    confirm.click();
    expect(onConfirm).toHaveBeenCalledTimes(1);
  });
});
