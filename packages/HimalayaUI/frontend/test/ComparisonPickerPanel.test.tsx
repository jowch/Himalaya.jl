import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, vi } from "vitest";
import { ComparisonPickerPanel } from "../src/components/ComparisonPickerPanel";
import * as api from "../src/api";
import { useAppState } from "../src/state";

const ROW = {
  sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
  indexing_exposure_id: 100,
  all_exposures: [{ id: 100, sample_id: 10, filename: "f1.dat", selected: true }],
};

beforeEach(() => {
  // Reset Zustand to a fresh draft so per-test seeding (or absence of) is
  // deterministic. The store reads `loadDraftFromSession` at module load,
  // so without an explicit reset prior tests' drafts leak.
  useAppState.getState().discardDraft();
  useAppState.getState().startNewDraft();

  vi.spyOn(api, "getPickerSamples").mockResolvedValue([ROW]);
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([]);
});
afterEach(() => { vi.restoreAllMocks(); });

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

test("inline shell does not render dialog role or focus trap", async () => {
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  expect(screen.queryByRole("dialog")).toBeNull();
});

test("toggle fires addMember immediately (immediate-commit semantics)", async () => {
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  // After immediate-commit, the draft should contain a member with exposure_id 100.
  const draft = useAppState.getState().activeDraft;
  expect(draft?.members.some((m) => m.exposure_id === 100)).toBe(true);
});

test("already-added rows render locked (alreadyAddedExposureIds gate)", async () => {
  // Pre-seed a draft with exposure 100 already a member.
  useAppState.getState().startNewDraft();
  // QueryClient mock so addMember's snapshot derivation doesn't blow up.
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  useAppState.getState().addMember(100, qc);
  wrap(<ComparisonPickerPanel experimentId={1} />);
  await screen.findByText("S1");
  // Re-clicking should not append a duplicate member (panel filters via
  // alreadyAddedExposureIds before calling addMember).
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  const members = useAppState.getState().activeDraft?.members ?? [];
  expect(members.filter((m) => m.exposure_id === 100).length).toBe(1);
});
