/**
 * Shell-level tests for ComparisonPicker.
 *
 * Body-level concerns (filter chips, recents, list rendering, locking) live
 * in ComparisonPickerBody.test.tsx. This file asserts only what the modal
 * shell owns: overlay visibility, Esc-to-close, outside-click, focus trap,
 * and the "Add N selected" → addMember integration.
 */
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, describe, it, expect, vi } from "vitest";
import userEvent from "@testing-library/user-event";
import { ComparisonPicker } from "../src/components/ComparisonPicker";
import { useAppState } from "../src/state";
import * as api from "../src/api";

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

function resetStore(): void {
  localStorage.clear();
  sessionStorage.clear();
  useAppState.setState({
    username: "alice",
    activeExperimentId: 1,
    activeDraft: {
      id: undefined,
      baseHash: undefined,
      title: "",
      description: "",
      members: [],
      forkedFromId: undefined,
      forkedAtHash: undefined,
    },
  });
}

beforeEach(() => {
  resetStore();
  // Minimal mocks so the body doesn't throw during shell tests.
  vi.spyOn(api, "getPickerSamples").mockResolvedValue([
    {
      sample: { id: 10, experiment_id: 1, name: "S1", display_name: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [{ id: 100, sample_id: 10, filename: "f1.dat", selected: true }],
    },
    {
      sample: { id: 20, experiment_id: 1, name: "S2", display_name: null, notes: null, tags: [] },
      indexing_exposure_id: 200,
      all_exposures: [{ id: 200, sample_id: 20, filename: "f2.dat", selected: true }],
    },
  ]);
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([]);
  vi.spyOn(api, "listExperiments").mockResolvedValue([
    { id: 1, name: "E", config: null },
  ]);
});

afterEach(() => { vi.restoreAllMocks(); });

describe("<ComparisonPicker> — shell", () => {
  it("not rendered when isOpen=false", () => {
    wrap(<ComparisonPicker isOpen={false} onClose={() => {}} experimentId={1} />);
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("Esc key fires onClose", () => {
    const onClose = vi.fn();
    wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
    fireEvent.keyDown(screen.getByRole("dialog"), { key: "Escape" });
    expect(onClose).toHaveBeenCalled();
  });

  it("clicking outside fires onClose", () => {
    const onClose = vi.fn();
    wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
    fireEvent.click(screen.getByTestId("comparison-picker-overlay"));
    expect(onClose).toHaveBeenCalled();
  });

  it("title is 'Add traces'", () => {
    wrap(<ComparisonPicker isOpen={true} onClose={() => {}} experimentId={1} />);
    expect(screen.getByText("Add traces")).toBeInTheDocument();
  });

  it("close button fires onClose", () => {
    const onClose = vi.fn();
    wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
    fireEvent.click(screen.getByTestId("comparison-picker-close"));
    expect(onClose).toHaveBeenCalled();
  });

  it("cancel button fires onClose", () => {
    const onClose = vi.fn();
    wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);
    fireEvent.click(screen.getByTestId("comparison-picker-cancel"));
    expect(onClose).toHaveBeenCalled();
  });

  it("'Add selected' button is disabled when nothing is picked", () => {
    wrap(<ComparisonPicker isOpen={true} onClose={() => {}} experimentId={1} />);
    const addBtn = screen.getByTestId("comparison-picker-add");
    expect((addBtn as HTMLButtonElement).disabled).toBe(true);
  });

  it("search input is focused on cold open — regression: PR #96 review", () => {
    // Override the beforeEach's mockResolvedValue with a never-resolving
    // promise so pickerQ.isLoading stays true. If a future refactor moves
    // the search input back inside the Skeleton boundary, inputRef.current
    // is null at effect time and this assertion breaks.
    vi.spyOn(api, "getPickerSamples").mockReturnValue(new Promise<never>(() => {}));
    wrap(<ComparisonPicker isOpen={true} onClose={() => {}} experimentId={1} />);
    expect(document.activeElement).toBe(
      screen.getByTestId("comparison-picker-search"),
    );
  });

  it("focus trap: Tab cycles within dialog", async () => {
    const user = userEvent.setup();
    wrap(<ComparisonPicker isOpen={true} onClose={() => {}} experimentId={1} />);
    const dialog = screen.getByRole("dialog");
    const focusable = dialog.querySelectorAll(
      'button:not([disabled]),input:not([disabled]),[tabindex]:not([tabindex="-1"])',
    );
    expect(focusable.length).toBeGreaterThan(1);
    const first = focusable[0] as HTMLElement;
    const last = focusable[focusable.length - 1] as HTMLElement;
    last.focus();
    expect(document.activeElement).toBe(last);
    await user.keyboard("{Tab}");
    expect(document.activeElement).toBe(first);
  });

  it("'Add N selected' fires addMember per pick on click", async () => {
    const user = userEvent.setup();
    const onClose = vi.fn();
    wrap(<ComparisonPicker isOpen={true} onClose={onClose} experimentId={1} />);

    // Wait for body rows to appear.
    const s1Checkbox = await screen.findByText("S1")
      .then((el) => el.closest("[data-testid='sample-picker-row']")!)
      .then((row) => row.querySelector("[data-testid='sample-picker-row-checkbox']") as HTMLInputElement);
    await user.click(s1Checkbox);

    // Button should now say "Add 1 selected".
    const addBtn = screen.getByTestId("comparison-picker-add");
    expect(addBtn).toHaveTextContent("Add 1 selected");

    await user.click(addBtn);

    // onClose was called.
    expect(onClose).toHaveBeenCalledTimes(1);

    // activeDraft should have one member with exposure_id 100.
    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    expect(draft!.members).toHaveLength(1);
    expect(draft!.members[0]!.exposure_id).toBe(100);
  });
});
