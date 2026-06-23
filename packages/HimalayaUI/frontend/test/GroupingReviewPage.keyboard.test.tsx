import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

// Keyboard nav on the grouping review: ↑/↓ move a roving cursor through the
// visible samples, Space toggles selection of the cursored sample, Shift+↑/↓
// jump to the prev/next flagged sample, and x dismisses the cursored flag.

const dismissMutate = vi.fn();

// One load, 3 samples: S1 (clean), S2 (split-flagged), S3 (clean).
const LOADS: Load[] = [
  { load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 1, name: "S1", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 2, name: "S2", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "split", split_at_index: 2, jump_from: 8, jump_to: 36 }, exposures: [] },
      { sample_id: 3, name: "S3", slot_index: 3, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
    ] },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: () => ({ mutate: vi.fn(), isPending: false }),
    useRenameSample: () => ({ mutate: vi.fn(), isPending: false }),
    useMoveExposure: () => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: () => ({ mutate: vi.fn(), isPending: false }),
    useDismissGroupingFlag: () => ({ mutate: dismissMutate, isPending: false }),
    useUndoDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
  };
});

vi.mock("../src/lib/toast", () => ({ showToast: vi.fn() }));

function wrap(node: ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

/** Show ALL samples (the default "Needs review" filter would hide the clean ones). */
function showAll() {
  fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
}

/** The sample-fold for a given name. */
function foldFor(name: string): HTMLElement {
  const heading = screen.getByText(name, { selector: "span" });
  return heading.closest('[data-testid="sample-fold"]') as HTMLElement;
}

beforeEach(() => { dismissMutate.mockClear(); });

describe("GroupingReviewPage keyboard nav", () => {
  it("ArrowDown moves the cursor onto the first sample, then the next", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(foldFor("S1")).toHaveAttribute("data-cursored", "true");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
    expect(foldFor("S1")).not.toHaveAttribute("data-cursored");
  });

  it("Space toggles selection of the cursored sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    fireEvent.keyDown(window, { key: "ArrowDown" }); // cursor → S1
    fireEvent.keyDown(window, { key: " " });
    expect(within(foldFor("S1")).getByRole("checkbox")).toBeChecked();
  });

  it("Shift+ArrowDown jumps the cursor to the next flagged sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    fireEvent.keyDown(window, { key: "ArrowDown" }); // cursor → S1 (clean)
    fireEvent.keyDown(window, { key: "ArrowDown", shiftKey: true }); // → S2 (flagged)
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("x dismisses the flag on the cursored sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    fireEvent.keyDown(window, { key: "ArrowDown", shiftKey: true }); // cursor → S2 (first flagged)
    fireEvent.keyDown(window, { key: "x" });
    expect(dismissMutate).toHaveBeenCalledTimes(1);
    expect(dismissMutate.mock.calls[0]![0]).toMatchObject({ sampleId: 2 });
  });

  it("x on a clean (unflagged) cursored sample does nothing", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    fireEvent.keyDown(window, { key: "ArrowDown" }); // cursor → S1 (clean)
    fireEvent.keyDown(window, { key: "x" });
    expect(dismissMutate).not.toHaveBeenCalled();
  });
});
