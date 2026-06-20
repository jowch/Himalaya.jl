import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const mergeMutate = vi.fn();
const LOADS: Load[] = [
  { load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 11, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 12, name: "C", slot_index: 3, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
    ] },
];
vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  // Hooks take only experimentId (settled in Task 9 — entity ids go in mutate input).
  return {
    ...actual, useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: (_experimentId: number) => ({ mutate: mergeMutate, isPending: false }),
    useRenameSample: (_experimentId: number) => ({ mutate: vi.fn() }),
    useMoveExposure: (_experimentId: number) => ({ mutate: vi.fn() }),
    useSplitSample: (_experimentId: number) => ({ mutate: vi.fn() }),
    useDismissGroupingFlag: (_experimentId: number) => ({ mutate: vi.fn() }),
  };
});
const wrap = (n: ReactNode) => render(<QueryClientProvider client={new QueryClient()}>{n}</QueryClientProvider>);
beforeEach(() => mergeMutate.mockClear());

describe("bulk merge", () => {
  it("merges all selected losers into the first-selected survivor (one mutate per loser)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select A/i }));   // survivor
    fireEvent.click(screen.getByRole("checkbox", { name: /select B/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select C/i }));
    // bulk-bar primary (scope the query to the bar so it doesn't collide with
    // the later modal confirm button, which also reads "Merge")
    const bar = screen.getByTestId("bulk-bar");
    fireEvent.click(within(bar).getByRole("button", { name: /^merge$/i }));
    fireEvent.click(screen.getByRole("button", { name: /^merge$/i }));        // confirm modal
    expect(mergeMutate).toHaveBeenCalledTimes(2);
    // assert first-arg inputs (the page may pass a 2nd options arg)
    const inputs = mergeMutate.mock.calls.map((c) => c[0]);
    expect(inputs).toContainEqual({ loserId: 11, survivorId: 10 });
    expect(inputs).toContainEqual({ loserId: 12, survivorId: 10 });
  });
});
