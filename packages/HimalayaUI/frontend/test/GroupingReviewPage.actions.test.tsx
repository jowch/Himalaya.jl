// test/GroupingReviewPage.actions.test.tsx
//
// Task 6.1 — GroupingReviewPage interaction-architecture migration.
// Mounts the full TestShell (keyboard layer + InteractionDock) so dock buttons
// AND keyboard-triggered actions flow through the REAL registry path.
// Arrow keys (↑/↓) fire on the scope container (data-interaction-scope).
// Space, x, Shift+Arrow fire on window (keyboard layer).
import { describe, it, expect, vi, beforeEach, type Mock } from "vitest";
import type React from "react";
import { render, screen, fireEvent, within, act } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/pages/GroupingReviewPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { useListCursor } from "../src/print/interaction/useListCursor";
import { DockStepper } from "../src/print/ui";
import { runCursorContract } from "./interaction/cursorContract";

// ── fixtures ──────────────────────────────────────────────────────────────────
// THREE SEPARATE LOADS — one sample each. attn filter keeps only loads that
// contain a flagged sample, so attn-filter → Load 2 only → flatSampleIds=[2]
// → cursor initializes on S2. showAll() exposes [S1,S2,S3]; cursor stays on S2
// (still in list, index 1) for ArrowUp→S1 / ArrowDown→S3 tests.
const IDS = [1, 2, 3];
const LOADS: Load[] = [
  {
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
    frame_count: 0, note: null,
    samples: [
      { sample_id: 1, name: "S1", slot_index: 1, grouping_source: "auto_position", name_source: "auto",
        merged_into_id: null, flag: null, exposures: [] },
    ],
  },
  {
    load_id: 2, load_index: 2, session_id: null, start_time: null, end_time: null,
    frame_count: 0, note: null,
    samples: [
      { sample_id: 2, name: "S2", slot_index: 2, grouping_source: "auto_position", name_source: "auto",
        merged_into_id: null, flag: { kind: "split", split_at_index: 2, jump_from: 8, jump_to: 36 }, exposures: [] },
    ],
  },
  {
    load_id: 3, load_index: 3, session_id: null, start_time: null, end_time: null,
    frame_count: 0, note: null,
    samples: [
      { sample_id: 3, name: "S3", slot_index: 3, grouping_source: "auto_position", name_source: "auto",
        merged_into_id: null, flag: null, exposures: [] },
    ],
  },
];

const splitMutate = vi.fn();
const dismissMutate = vi.fn();

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: () => ({ mutate: vi.fn(), isPending: false }),
    useRenameSample: () => ({ mutate: vi.fn(), isPending: false }),
    useMoveExposure: () => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: () => ({ mutate: splitMutate, isPending: false }),
    useDismissGroupingFlag: () => ({ mutate: dismissMutate, isPending: false }),
    useUndoDismissGroupingFlag: () => ({ mutate: vi.fn(), isPending: false }),
  };
});

vi.mock("../src/lib/toast", () => ({ showToast: vi.fn() }));

// ── TestShell ─────────────────────────────────────────────────────────────────
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return (
    <>
      {children}
      <InteractionDock />
    </>
  );
}

function renderGrouping(props?: { onBack?: () => void; onConfirm?: () => void }) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <TestShell>
        <GroupingReviewPage experimentId={7} onBack={props?.onBack ?? vi.fn()} onConfirm={props?.onConfirm} />
      </TestShell>
    </QueryClientProvider>,
  );
}

function getScope(): HTMLElement {
  const scope = document.querySelector("[data-interaction-scope]") as HTMLElement | null;
  if (!scope) throw new Error("No [data-interaction-scope] found");
  return scope;
}

function foldFor(name: string): HTMLElement {
  const heading = screen.getByText(name, { selector: "span" });
  return heading.closest('[data-testid="sample-fold"]') as HTMLElement;
}

function showAll() {
  fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
}

beforeEach(() => {
  dismissMutate.mockClear();
  splitMutate.mockClear();
  useInteraction.getState().clearPage();
});

// ── tests ─────────────────────────────────────────────────────────────────────
describe("GroupingReviewPage interaction", () => {
  it("cursor initializes on the first visible sample (S2 in attn filter)", () => {
    renderGrouping();
    // attn filter: only S2 (flagged) is visible → cursor lands on S2
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("ArrowDown on scope advances cursor (S2→S3 after showAll)", () => {
    renderGrouping();
    showAll(); // cursor still on S2; now S1,S2,S3 all visible (S2 at index 1)
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowDown" }); });
    expect(foldFor("S3")).toHaveAttribute("data-cursored", "true");
    expect(foldFor("S2")).not.toHaveAttribute("data-cursored");
  });

  it("ArrowUp on scope retreats cursor (S2→S1 after showAll)", () => {
    renderGrouping();
    showAll(); // cursor on S2 (index 1)
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); });
    expect(foldFor("S1")).toHaveAttribute("data-cursored", "true");
  });

  it("Shift+ArrowDown does NOT move the cursor (routes to nextFlagged action instead)", () => {
    renderGrouping();
    showAll();
    // cursor on S2; Shift+ArrowDown fires on scope but !e.shiftKey guard passes it
    // to the keyboard layer as nextFlagged — it does not become a cursor.moveBy call
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowDown", shiftKey: true }); });
    // cursor should still be on S2 (nextFlagged from S2 finds nothing after)
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("Space (keyboard layer) toggles the page's ordered selection on the cursored sample", () => {
    renderGrouping();
    // cursor on S2 in attn filter — no need for showAll
    act(() => { fireEvent.keyDown(window, { key: " " }); });
    expect(within(foldFor("S2")).getByRole("checkbox")).toBeChecked();
    // Space again deselects
    act(() => { fireEvent.keyDown(window, { key: " " }); });
    expect(within(foldFor("S2")).getByRole("checkbox")).not.toBeChecked();
  });

  it("selection order is preserved (first selected = survivor)", () => {
    renderGrouping();
    showAll();
    // S2 is at cursor; select it first
    act(() => { fireEvent.keyDown(window, { key: " " }); });
    // Move to S1 and select second
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); });
    act(() => { fireEvent.keyDown(window, { key: " " }); });
    // dock shows merge (2 selected, ≥2 threshold met)
    const mergeBtn = screen.getByTestId("dock-action-merge");
    expect(mergeBtn).toBeEnabled();
    // clicking merge opens bulk-merge-confirm with S2 as survivor
    fireEvent.click(mergeBtn);
    expect(screen.getByTestId("bulk-merge-confirm")).toBeInTheDocument();
    // Confirm modal text names S2 as the survivor
    expect(within(screen.getByTestId("bulk-merge-confirm")).getByText(/S2/)).toBeInTheDocument();
  });

  it("merge dock button hidden when no samples selected (mode: selection)", () => {
    renderGrouping();
    showAll();
    // Nothing selected → merge button not rendered (mode:selection hides it)
    expect(screen.queryByTestId("dock-action-merge")).toBeNull();
  });

  it("merge dock button appears and is enabled at ≥2 selected", () => {
    renderGrouping();
    showAll();
    // Select S2 via checkbox, then S1 via checkbox
    fireEvent.click(within(foldFor("S2")).getByRole("checkbox"));
    fireEvent.click(within(foldFor("S1")).getByRole("checkbox"));
    const mergeBtn = screen.getByTestId("dock-action-merge");
    expect(mergeBtn).toBeEnabled();
    fireEvent.click(mergeBtn);
    expect(screen.getByTestId("bulk-merge-confirm")).toBeInTheDocument();
  });

  it("Shift+ArrowDown (keyboard layer) jumps cursor to the next flagged sample", () => {
    renderGrouping();
    showAll();
    // cursor on S2; move to S1 first
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); });
    expect(foldFor("S1")).toHaveAttribute("data-cursored", "true");
    // nextFlagged: from S1 → S2
    act(() => { fireEvent.keyDown(window, { key: "ArrowDown", shiftKey: true }); });
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("Shift+ArrowUp (keyboard layer) jumps cursor to the prev flagged sample", () => {
    renderGrouping();
    showAll();
    // cursor on S2; move to S3 first
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowDown" }); });
    expect(foldFor("S3")).toHaveAttribute("data-cursored", "true");
    // prevFlagged: from S3 → S2
    act(() => { fireEvent.keyDown(window, { key: "ArrowUp", shiftKey: true }); });
    expect(foldFor("S2")).toHaveAttribute("data-cursored", "true");
  });

  it("x (keyboard layer) dismisses the flag on the cursored flagged sample", () => {
    renderGrouping();
    // cursor starts on S2 (flagged) in attn filter
    act(() => { fireEvent.keyDown(window, { key: "x" }); });
    expect(dismissMutate).toHaveBeenCalledTimes(1);
    expect(dismissMutate.mock.calls[0]![0]).toMatchObject({ sampleId: 2 });
  });

  it("x on a clean (unflagged) cursored sample does not call dismissMutate", () => {
    renderGrouping();
    showAll();
    act(() => { fireEvent.keyDown(getScope(), { key: "ArrowUp" }); }); // cursor → S1 (clean)
    act(() => { fireEvent.keyDown(window, { key: "x" }); });
    expect(dismissMutate).not.toHaveBeenCalled();
  });

  it("x typed in a rename input (isTyping guard) does NOT dismiss a flag", () => {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    function Shell({ children }: { children: React.ReactNode }): JSX.Element {
      useKeyboardLayer();
      return (
        <>
          {children}
          <input data-testid="rename-input" defaultValue="" />
          <InteractionDock />
        </>
      );
    }
    render(
      <QueryClientProvider client={qc}>
        <Shell>
          <GroupingReviewPage experimentId={7} onBack={vi.fn()} />
        </Shell>
      </QueryClientProvider>,
    );
    const input = screen.getByTestId("rename-input");
    act(() => { input.focus(); });
    act(() => { fireEvent.keyDown(input, { key: "x" }); });
    expect(dismissMutate).not.toHaveBeenCalled();
  });

  it("confirm dock button is enabled when not scanning", () => {
    renderGrouping();
    expect(screen.getByTestId("dock-action-confirm")).toBeEnabled();
  });

  it("confirm dock button fires onConfirm when provided", () => {
    const onConfirm = vi.fn();
    renderGrouping({ onConfirm });
    fireEvent.click(screen.getByTestId("dock-action-confirm"));
    expect(onConfirm).toHaveBeenCalledTimes(1);
  });

  it("back up-link calls onBack", () => {
    const onBack = vi.fn();
    renderGrouping({ onBack });
    fireEvent.click(screen.getByTestId("dock-up-link"));
    expect(onBack).toHaveBeenCalledTimes(1);
  });

  it("select dock button toggles selection via click", () => {
    renderGrouping();
    // cursor on S2 (attn filter); select dock button is enabled
    const selectBtn = screen.getByTestId("dock-action-select");
    expect(selectBtn).toBeEnabled();
    fireEvent.click(selectBtn);
    expect(within(foldFor("S2")).getByRole("checkbox")).toBeChecked();
  });
});

// ── cursor contract (headless probe) ─────────────────────────────────────────
// GroupingReviewPage uses a HEADLESS cursor (no rowProps on LoadFold rows).
// The contract runs on a standalone Probe that DOES spread rowProps — this
// validates the cursor's own contract, independent of the page's headless usage.
runCursorContract("Grouping sample cursor", () => {
  const onActivate = vi.fn() as Mock;
  return {
    ui: (capture) => {
      function Probe(): JSX.Element {
        const curs = useListCursor({
          ids: IDS,
          onActivate,
          stepperLabel: "Sample",
          stepperTestIdBase: "sample",
          axis: "vertical",
        });
        capture({ cursor: curs, ids: IDS, onActivate });
        return (
          <div>
            {IDS.map((id) => (
              <div key={id} {...curs.rowProps(id)}>
                sample {id}
              </div>
            ))}
            <DockStepper {...curs.stepperProps()} />
          </div>
        );
      }
      return <Probe />;
    },
  };
});
