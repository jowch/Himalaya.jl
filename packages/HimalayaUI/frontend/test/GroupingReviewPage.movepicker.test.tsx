/**
 * GroupingReviewPage — Move picker (same-load sibling Menu).
 *
 * FIX (b): activating Move on an exposure must open a Menu listing every
 * OTHER sample in the same load (not the current owner). Selecting one fires
 * moveMutate({ exposureId, sampleId: <destination> }) where destination !=
 * current owner. A self-move must NEVER be issued.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const moveMutate = vi.fn();

// Two loads:
//   Load 1: samples S10 (owns exp E100) and S11 (sibling in same load)
//   Load 2: sample S20 only (no siblings — Move picker should show no options)
const LOADS: Load[] = [
  {
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
    frame_count: 0, note: null,
    samples: [
      {
        sample_id: 10, name: "Alpha", slot_index: 1, grouping_source: "auto_position",
        name_source: "auto", merged_into_id: null,
        flag: { kind: "split", split_at_index: 1, jump_from: 1.0, jump_to: 5.0 },
        exposures: [
          { id: 100, filename: "img_001.tif", horizontal_position: 1.0, timestamp: "10:00:00" },
        ],
      },
      {
        sample_id: 11, name: "Beta", slot_index: 2, grouping_source: "auto_position",
        name_source: "auto", merged_into_id: null, flag: null,
        exposures: [],
      },
    ],
  },
  {
    load_id: 2, load_index: 2, session_id: null, start_time: null, end_time: null,
    frame_count: 0, note: null,
    samples: [
      {
        sample_id: 20, name: "Gamma", slot_index: 1, grouping_source: "auto_position",
        name_source: "auto", merged_into_id: null, flag: null,
        exposures: [
          { id: 200, filename: "img_002.tif", horizontal_position: 2.0, timestamp: "11:00:00" },
        ],
      },
    ],
  },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: () => ({ mutate: vi.fn(), isPending: false }),
    useRenameSample: () => ({ mutate: vi.fn(), isPending: false }),
    useMoveExposure: () => ({ mutate: moveMutate, isPending: false }),
    useSplitSample: () => ({ mutate: vi.fn(), isPending: false }),
    useDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
    useUndoDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
  };
});

// The open folds render exposure thumbnails through DetectorImage, which
// touches fetch / createImageBitmap (absent in JSDOM).
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

vi.mock("../src/lib/toast", () => ({
  showToast: vi.fn(),
}));

function wrap(node: ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

beforeEach(() => {
  moveMutate.mockClear();
});

describe("GroupingReviewPage — Move picker", () => {
  it("activating Move on an exposure opens a Menu of same-load siblings (NOT the current owner)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);

    // Switch to "All loads" so Load 1 (split flag) and Load 2 are both visible.
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));

    // The Move button on exposure E100 (owned by S10=Alpha, sibling S11=Beta).
    const moveBtn = screen.getByRole("button", { name: /move to another sample/i });
    fireEvent.click(moveBtn);

    // The picker menu must be open.
    const menu = screen.getByRole("menu");
    expect(menu).toBeInTheDocument();

    // The menu must list Beta (sibling S11) but NOT Alpha (current owner S10).
    expect(screen.getByRole("menuitem", { name: /Beta/i })).toBeInTheDocument();
    expect(screen.queryByRole("menuitem", { name: /Alpha/i })).toBeNull();
  });

  it("selecting a sibling from the picker fires moveMutate with the destination sampleId (not the current owner)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);

    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));

    const moveBtn = screen.getByRole("button", { name: /move to another sample/i });
    fireEvent.click(moveBtn);

    // Select Beta (sample_id=11).
    fireEvent.click(screen.getByRole("menuitem", { name: /Beta/i }));

    expect(moveMutate).toHaveBeenCalledTimes(1);
    const call = moveMutate.mock.calls[0]![0];
    expect(call).toMatchObject({ exposureId: 100, sampleId: 11 });
    // Crucially: destination must differ from current owner (10).
    expect(call.sampleId).not.toBe(10);
  });

  it("after selecting a destination the picker menu is dismissed", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);

    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));

    const moveBtn = screen.getByRole("button", { name: /move to another sample/i });
    fireEvent.click(moveBtn);

    fireEvent.click(screen.getByRole("menuitem", { name: /Beta/i }));

    // Menu must be gone after selection.
    expect(screen.queryByRole("menu")).toBeNull();
  });
});
