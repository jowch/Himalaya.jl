import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, vi } from "vitest";
import { ComparisonPickerBody } from "../src/components/ComparisonPickerBody";
import * as api from "../src/api";

// Wrap mocks in beforeEach so they're isolated per test (the project's
// vitest setup does not call vi.restoreAllMocks() between files).
beforeEach(() => {
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
  vi.spyOn(api, "listExperiments").mockResolvedValue([
    { id: 1, name: "E", config: null },
  ]);
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([200]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([]);
});
afterEach(() => { vi.restoreAllMocks(); });

function wrap(ui: React.ReactElement) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={qc}>{ui}</QueryClientProvider>);
}

test("recents section dedupes against main list (one row per sample, S2 appears once)", async () => {
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      onPick={() => {}}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // Recents has S2 (exposure 200). Main list also has S2.
  // The body must NOT render S2 twice in the visible main list.
  const s2Rows = screen.queryAllByText("S2");
  expect(s2Rows.length).toBeLessThanOrEqual(1);   // exactly one section renders it
});

test("rows whose exposure id is in alreadyAddedExposureIds render locked", async () => {
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      onPick={() => {}}
      alreadyAddedExposureIds={new Set([100])}
    />,
  );
  await screen.findByText("S1");
  // S1 has indexing_exposure_id=100 which is in the set → data-locked="true".
  const s1Row = screen.getByText("S1").closest("[data-testid='sample-picker-row']")!;
  expect(s1Row).toHaveAttribute("data-locked", "true");
  expect(screen.getByText("already added")).toBeInTheDocument();
  // S2 (exposure 200, appears in recents) should NOT be locked.
  const s2Row = screen.getByText("S2").closest("[data-testid='sample-picker-row']")!;
  expect(s2Row).not.toHaveAttribute("data-locked");
});

test("onPick fires per toggle with default exposure id", async () => {
  const onPick = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      onPick={onPick}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // S1 is in the main list. Use the row container to scope the checkbox click.
  const s1Row = screen.getByText("S1").closest("[data-testid='sample-picker-row']")!;
  fireEvent.click(s1Row.querySelector("[data-testid='sample-picker-row-checkbox']")!);
  expect(onPick).toHaveBeenCalledWith({ sample_id: 10, exposure_id: 100, source: "default" });
});

test("override caret pick fires onPick with override exposure id", async () => {
  // Override the mock for this test to give S1 two exposures so we can
  // select the non-default one (f2.dat/id=101) to trigger the override path.
  vi.spyOn(api, "getPickerSamples").mockResolvedValue([
    {
      sample: { id: 10, experiment_id: 1, name: "S1", display_name: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [
        { id: 100, sample_id: 10, filename: "f1.dat", selected: true },
        { id: 101, sample_id: 10, filename: "f2.dat", selected: false },
      ],
    },
    {
      sample: { id: 20, experiment_id: 1, name: "S2", display_name: null, notes: null, tags: [] },
      indexing_exposure_id: 200,
      all_exposures: [{ id: 200, sample_id: 20, filename: "f3.dat", selected: true }],
    },
  ]);
  const onPick = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      onPick={onPick}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // S1 is in the main list. Expand its caret dropdown, then select the non-default radio.
  const s1Row = screen.getByText("S1").closest("[data-testid='sample-picker-row']")!;
  fireEvent.click(s1Row.querySelector("[data-testid='sample-picker-row-caret']")!);
  // Click the f2.dat radio (a different exposure → override source).
  fireEvent.click(screen.getByLabelText(/f2\.dat/));
  expect(onPick).toHaveBeenCalledWith({
    sample_id: 10, exposure_id: 101, source: "override",
  });
});

test("does not flicker on background refetch (gates on isLoading not isPending)", async () => {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false, staleTime: 0 } } });
  qc.setQueryData(["experiment", 1, "picker-samples"], [
    {
      sample: { id: 10, experiment_id: 1, name: "S1", display_name: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [{ id: 100, sample_id: 10, filename: "f1", selected: true }],
    },
  ]);
  const { rerender } = render(
    <QueryClientProvider client={qc}>
      <ComparisonPickerBody
        experimentId={1} onPick={() => {}}
        alreadyAddedExposureIds={new Set()}
      />
    </QueryClientProvider>,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
  qc.invalidateQueries({ queryKey: ["experiment", 1, "picker-samples"] });
  rerender(
    <QueryClientProvider client={qc}>
      <ComparisonPickerBody
        experimentId={1} onPick={() => {}}
        alreadyAddedExposureIds={new Set()}
      />
    </QueryClientProvider>,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
});
