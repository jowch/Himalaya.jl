import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

// Regression: agbe (2 exposures) had a 1-BASED split_at_index of 2; the page
// treated it as 0-based and sliced exposures.slice(2) === [] → the backend
// rejected the split with "exposure_ids must not be empty". The fix maps
// split_at_index - 1 (clamped to [1, len-1]) so the split-off tail is the
// 0-based index-1 exposure onward — never empty.

const splitMutate = vi.fn();

function loadsWith(splitAtIndex: number): Load[] {
  return [
    { load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
      samples: [{ sample_id: 30, name: "agbe (S01P01)", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "split", split_at_index: splitAtIndex, jump_from: 54.4, jump_to: 43.0 },
        exposures: [
          { id: 100, filename: "agbe.tif", horizontal_position: 54.4, timestamp: "10:00" },
          { id: 101, filename: "emptycap.tif", horizontal_position: 43.0, timestamp: "10:01" },
        ] }] },
  ];
}

let LOADS: Load[] = loadsWith(2);

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: () => ({ mutate: vi.fn(), isPending: false }),
    useRenameSample: () => ({ mutate: vi.fn(), isPending: false }),
    useMoveExposure: () => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: () => ({ mutate: splitMutate, isPending: false }),
    useDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
    useUndoDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
  };
});

vi.mock("../src/lib/toast", () => ({ showToast: vi.fn() }));

// DetectorImage touches fetch / createImageBitmap (absent in JSDOM); the open
// split fold now renders exposure thumbnails through it.
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

function wrap(node: ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

beforeEach(() => { splitMutate.mockClear(); LOADS = loadsWith(2); });

describe("GroupingReviewPage split", () => {
  it("a 1-based split_at_index=2 on a 2-exposure sample splits off the 2nd exposure (never empty)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: "Split..." }));
    expect(splitMutate).toHaveBeenCalledTimes(1);
    const arg = splitMutate.mock.calls[0]![0] as { sampleId: number; exposureIds: number[] };
    expect(arg.sampleId).toBe(30);
    expect(arg.exposureIds).toEqual([101]); // the split-off tail, non-empty
  });

  it("clamps a degenerate split_at_index past the end so exposure_ids is never empty", () => {
    LOADS = loadsWith(3); // 1-based 3 on a 2-exposure sample → would slice past the end
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: "Split..." }));
    const arg = splitMutate.mock.calls[0]![0] as { exposureIds: number[] };
    expect(arg.exposureIds.length).toBeGreaterThan(0);
    expect(arg.exposureIds).toEqual([101]);
  });
});
