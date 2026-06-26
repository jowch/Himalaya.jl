// test/GroupingReviewPage.keyboard.test.tsx
//
// Keyboard nav on the grouping review (migrated to interaction-arch model):
// ↑/↓ on the scope container move a roving cursor through visible samples.
// Space (keyboard layer) toggles selection. Shift+↑/↓ (keyboard layer) jump
// to the prev/next flagged. x (keyboard layer) dismisses the cursored flag.
//
// CURSOR INIT: useListCursor initializes to ids[0]. With the "attn" filter
// (default), only S2 (flagged) is visible → cursor starts on S2.
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, within, act } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/pages/GroupingReviewPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";

const dismissMutate = vi.fn();

// THREE SEPARATE LOADS — one sample each. The attn filter keeps only loads that
// contain a flagged sample, so: attn-filter → Load 2 only → flatSampleIds=[2]
// → cursor initializes on S2. showAll() exposes [S1,S2,S3]; cursor stays on S2
// (still in list at index 1) enabling ArrowUp→S1 / ArrowDown→S3.
const LOADS: Load[] = [
  {
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 1, name: "S1", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
    ],
  },
  {
    load_id: 2, load_index: 2, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 2, name: "S2", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "split", split_at_index: 2, jump_from: 8, jump_to: 36 }, exposures: [] },
    ],
  },
  {
    load_id: 3, load_index: 3, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 3, name: "S3", slot_index: 3, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
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
    useMoveExposure: () => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: () => ({ mutate: vi.fn(), isPending: false }),
    useDismissGroupingFlag: () => ({ mutate: dismissMutate, isPending: false }),
    useUndoDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
  };
});

vi.mock("../src/lib/toast", () => ({ showToast: vi.fn() }));

// TestShell: mounts the keyboard layer (handles Space/x/Shift+Arrow) +
// InteractionDock (renders dock buttons for dock:true actions).
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return <>{children}<InteractionDock /></>;
}

function wrap(node: React.ReactNode) {
  return render(
    <QueryClientProvider client={new QueryClient()}><TestShell>{node}</TestShell></QueryClientProvider>,
  );
}

/** Show ALL samples (the default "Needs review" filter hides the clean ones). */
function showAll() {
  fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
}

/** The sample-fold for a given name. */
function foldFor(name: string): HTMLElement {
  const heading = screen.getByText(name, { selector: "span" });
  return heading.closest('[data-testid="sample-fold"]') as HTMLElement;
}

/** The scope container that handles Arrow key cursor navigation. */
function getScope(): HTMLElement {
  return document.querySelector("[data-interaction-scope]") as HTMLElement;
}

beforeEach(() => { dismissMutate.mockClear(); });

describe("GroupingReviewPage keyboard nav", () => {
  it("cursor initializes on the first visible sample (S2 in attn filter)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // attn filter shows only S2 (the only flagged sample)
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("ArrowDown on scope moves cursor from S2 to S3 (after showAll)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll(); // cursor stays on S2 (it's still in the id list); flatSamples = [S1,S2,S3]
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowDown" }); });
    // S2 is at index 1 → ArrowDown → S3 (index 2)
    expect(foldFor("S3")).toHaveAttribute("data-cursored", "true");
    expect(foldFor("S2")).not.toHaveAttribute("data-cursored");
  });

  it("ArrowUp on scope moves cursor from S2 to S1 (after showAll)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll(); // cursor on S2 (index 1 in [S1,S2,S3])
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); });
    expect(foldFor("S1")).toHaveAttribute("data-cursored", "true");
  });

  it("Space (keyboard layer) toggles selection of the cursored sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // cursor on S2 in attn filter — no need for showAll
    act(() => { fireEvent.keyDown(window, { key: " " }); });
    expect(within(foldFor("S2")).getByRole("checkbox")).toBeChecked();
  });

  it("Shift+ArrowDown (keyboard layer) jumps cursor to next flagged sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); }); // cursor → S1 (clean)
    act(() => { fireEvent.keyDown(window, { key: "ArrowDown", shiftKey: true }); }); // nextFlagged → S2
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("Shift+ArrowUp (keyboard layer) jumps cursor to prev flagged sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowDown" }); }); // cursor → S3 (clean)
    act(() => { fireEvent.keyDown(window, { key: "ArrowUp", shiftKey: true }); }); // prevFlagged → S2
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("x (keyboard layer) dismisses the flag on the cursored flagged sample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // cursor starts on S2 (flagged) in attn filter
    act(() => { fireEvent.keyDown(window, { key: "x" }); });
    expect(dismissMutate).toHaveBeenCalledTimes(1);
    expect(dismissMutate.mock.calls[0]![0]).toMatchObject({ sampleId: 2 });
  });

  it("x on a clean (unflagged) cursored sample does nothing", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    showAll();
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); }); // cursor → S1 (clean)
    act(() => { fireEvent.keyDown(window, { key: "x" }); });
    expect(dismissMutate).not.toHaveBeenCalled();
  });
});
