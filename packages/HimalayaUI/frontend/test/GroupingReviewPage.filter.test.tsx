import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const LOADS: Load[] = [
  { load_id: 1, load_index: 1, session_id: null, start_time: "10:02", end_time: "10:38", frame_count: 0, note: null,
    samples: [{ sample_id: 10, name: "HA85 (S01P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] }] },
  { load_id: 2, load_index: 2, session_id: null, start_time: "13:10", end_time: "13:51", frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "HA85 (S02P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA85 (S01P15)" }, exposures: [] }] },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return { ...actual, useLoads: () => ({ data: LOADS, isLoading: false }) };
});

function wrap(node: ReactNode) {
  const qc = new QueryClient();
  return render(<QueryClientProvider client={qc}>{node}</QueryClientProvider>);
}

describe("GroupingReviewPage filter + persistent selection", () => {
  it("defaults to 'Needs review' (only the flagged load 2 shows)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.queryByText("Load 1")).toBeNull();
    expect(screen.getByText("Load 2")).toBeInTheDocument();
  });

  it("'All loads' reveals every load", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    expect(screen.getByText("Load 1")).toBeInTheDocument();
    expect(screen.getByText("Load 2")).toBeInTheDocument();
  });

  it("selection PERSISTS across filter changes (the cross-load merge gesture)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // show all, select sample 10 in load 1
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select HA85 \(S01P15\)/i }));
    // search to sample 20's load only -- sample 10 leaves the view
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "S02" } });
    // the footer selection readout still reflects the 1 retained selection
    expect(screen.getByTestId("grouping-selection-count")).toHaveTextContent("1");
  });
});
