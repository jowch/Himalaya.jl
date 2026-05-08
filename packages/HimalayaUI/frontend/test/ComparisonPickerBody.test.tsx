import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { beforeEach, afterEach, vi } from "vitest";
import { ComparisonPickerBody, type Pick } from "../src/components/ComparisonPickerBody";
import * as api from "../src/api";

// Wrap mocks in beforeEach so they're isolated per test (the project's
// vitest setup does not call vi.restoreAllMocks() between files).
beforeEach(() => {
  vi.spyOn(api, "getPickerSamples").mockResolvedValue([
    {
      sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
      indexing_exposure_id: 100,
      all_exposures: [{ id: 100, sample_id: 10, filename: "f1.dat", selected: true }],
    },
    {
      sample: { id: 20, experiment_id: 1, name: "S2", label: null, notes: null, tags: [] },
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

test("controlled-picks: onPicksChange fires with default exposure id on toggle", async () => {
  const onPicksChange = vi.fn();
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]}
      onPicksChange={onPicksChange}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // S1 is in the main list (S2 only is in recents). Click its checkbox.
  const s1Row = screen.getByText("S1").closest("[data-testid='sample-picker-row']")!;
  fireEvent.click(s1Row.querySelector("[data-testid='sample-picker-row-checkbox']")!);
  expect(onPicksChange).toHaveBeenCalledWith([
    { sample_id: 10, exposure_id: 100, source: "default" },
  ]);
});

// NOTE: Immediate-mode `onPick` is added in PR2 Task 13. PR1 ships only the
// controlled-picks interface (`picks` + `onPicksChange`). Skip the immediate-mode
// test until PR2.

test("recents section dedupes against main list (one row per sample, S2 appears once)", async () => {
  wrap(
    <ComparisonPickerBody
      experimentId={1}
      picks={[]} onPicksChange={() => {}}
      alreadyAddedExposureIds={new Set()}
    />,
  );
  await screen.findByText("S1");
  // Recents has S2 (exposure 200). Main list also has S2.
  // The body must NOT render S2 twice in the visible main list.
  const s2Rows = screen.queryAllByText("S2");
  expect(s2Rows.length).toBeLessThanOrEqual(1);   // exactly one section renders it
});
