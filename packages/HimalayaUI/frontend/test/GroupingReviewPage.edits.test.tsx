import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/pages/GroupingReviewPage";

const mergeMutate = vi.fn();
const renameMutate = vi.fn();
const dismissMutate = vi.fn();
const undoDismissMutate = vi.fn();

const LOADS: Load[] = [
  { load_id: 2, load_index: 2, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "HA85 (S02P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA85 (S01P15)" }, exposures: [] }] },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  // Hooks now take only experimentId (settled in Task 9 — entity ids go in mutate input).
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: (_experimentId: number) => ({ mutate: mergeMutate, isPending: false }),
    useRenameSample: (_experimentId: number) => ({ mutate: renameMutate, isPending: false }),
    useMoveExposure: (_experimentId: number) => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: (_experimentId: number) => ({ mutate: vi.fn(), isPending: false }),
    useDismissGroupingFlag: (_experimentId: number) => ({ mutate: dismissMutate, isPending: false }),
    useUndoDismissGroupingFlag: (_experimentId: number) => ({ mutate: undoDismissMutate, isPending: false }),
  };
});

// Capture calls to showToast so we can exercise the Undo action slot.
const toastSpy = vi.fn();
vi.mock("../src/lib/toast", () => ({
  showToast: (...args: unknown[]) => toastSpy(...args),
}));

function wrap(node: ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

beforeEach(() => {
  mergeMutate.mockClear();
  renameMutate.mockClear();
  dismissMutate.mockClear();
  undoDismissMutate.mockClear();
  toastSpy.mockClear();
});

describe("GroupingReviewPage edits", () => {
  it("inline rename: activating rename and committing via Enter fires useRenameSample with trimmed name", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // Click the Rename button on the sample row
    fireEvent.click(screen.getByRole("button", { name: /rename sample/i }));
    // The inline input should appear
    const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
    fireEvent.change(input, { target: { value: "  HA85 Renamed  " } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(renameMutate).toHaveBeenCalledTimes(1);
    expect(renameMutate.mock.calls[0]![0]).toEqual({ sampleId: 20, name: "HA85 Renamed" });
    // Toast shown with undo
    expect(toastSpy).toHaveBeenCalledTimes(1);
  });

  it("inline rename: pressing Escape does NOT fire useRenameSample", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /rename sample/i }));
    const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
    fireEvent.change(input, { target: { value: "Cancelled" } });
    fireEvent.keyDown(input, { key: "Escape" });
    expect(renameMutate).not.toHaveBeenCalled();
  });

  it("clicking the merge-prompt's Merge opens a confirm, then fires useMergeSamples with {loserId,survivorId}", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /merge into that sample/i }));
    // confirm modal
    fireEvent.click(screen.getByRole("button", { name: /^merge$/i }));
    // assert on the FIRST arg only (the page may or may not pass a 2nd options arg)
    expect(mergeMutate).toHaveBeenCalledTimes(1);
    expect(mergeMutate.mock.calls[0]![0]).toEqual({ loserId: 20, survivorId: 10 });
  });

  it("'Keep separate' (dismiss) fires useDismissGroupingFlag and shows a toast with an Undo action that calls useUndoDismissGroupingFlag", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);

    // The sample has a merge flag, so the merge-prompt is open. Click "Keep separate".
    fireEvent.click(screen.getByRole("button", { name: /keep separate/i }));

    // dismiss mutator must have fired
    expect(dismissMutate).toHaveBeenCalledTimes(1);
    expect(dismissMutate.mock.calls[0]![0]).toMatchObject({ sampleId: 20 });

    // toast must have been shown
    expect(toastSpy).toHaveBeenCalledTimes(1);

    // The toast's action slot must trigger the undo mutator
    const toastArgs = toastSpy.mock.calls[0]!;
    const action = toastArgs[2] as { label: string; onClick: () => void } | undefined;
    expect(action).toBeDefined();
    expect(action!.label).toMatch(/undo/i);

    // Invoke the Undo action — undoDismissMutate must fire with sampleId
    action!.onClick();
    expect(undoDismissMutate).toHaveBeenCalledTimes(1);
    expect(undoDismissMutate.mock.calls[0]![0]).toMatchObject({ sampleId: 20 });
  });
});
