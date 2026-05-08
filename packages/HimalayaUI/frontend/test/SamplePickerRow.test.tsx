import { render, screen, fireEvent } from "@testing-library/react";
import { vi } from "vitest";
import { SamplePickerRow } from "../src/components/SamplePickerRow";
import type { PickerSampleRow } from "../src/api";

const baseRow: PickerSampleRow = {
  sample: { id: 10, experiment_id: 1, name: "S1", label: null, notes: null, tags: [] },
  indexing_exposure_id: 100,
  all_exposures: [
    { id: 100, sample_id: 10, filename: "f1.dat", selected: true },
    { id: 101, sample_id: 10, filename: "f2.dat", selected: false },
  ],
};

test("renders sample name as primary, no filename in primary slot", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByText("S1")).toBeInTheDocument();
  expect(screen.queryByText("f1.dat")).toBeNull();   // hidden until caret expanded
});

test("clicking row body toggles the checkbox — regression: PR #96 review (issue #53 carryover)", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={handle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  // Click the sample-name span — bubbles up to the outer div's onClick.
  fireEvent.click(screen.getByText("S1"));
  expect(handle).toHaveBeenCalledWith(true);
});

test("clicking the caret button does NOT toggle the checkbox", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={handle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  // Caret is a <button> — handleRowClick bails on closest("button").
  fireEvent.click(screen.getByTestId("sample-picker-row-caret"));
  expect(handle).not.toHaveBeenCalled();
});

test("clicking the checkbox does NOT double-fire via row bubble", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={handle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  // Checkbox is an HTMLInputElement — handleRowClick bails on it.
  // The change handler still fires for the actual checkbox toggle (once).
  fireEvent.click(screen.getByTestId("sample-picker-row-checkbox"));
  expect(handle).toHaveBeenCalledTimes(1);
});

test("disabled and alreadyAdded rows do not toggle on body click", () => {
  // Disabled (zero-exposure) row.
  const disabledHandle = vi.fn();
  const disabled: PickerSampleRow = {
    ...baseRow, indexing_exposure_id: null, all_exposures: [],
  };
  const { rerender } = render(
    <SamplePickerRow
      row={disabled} checked={false} onCheckedChange={disabledHandle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  fireEvent.click(screen.getByText("S1"));
  expect(disabledHandle).not.toHaveBeenCalled();

  // alreadyAdded row.
  const lockedHandle = vi.fn();
  rerender(
    <SamplePickerRow
      row={baseRow} checked={true} onCheckedChange={lockedHandle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
      alreadyAdded={true}
    />,
  );
  fireEvent.click(screen.getByText("S1"));
  expect(lockedHandle).not.toHaveBeenCalled();
});

test("caret toggles override list visibility", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  const caret = screen.getByTestId("sample-picker-row-caret");
  fireEvent.click(caret);
  expect(screen.getByText("f1.dat")).toBeInTheDocument();
  expect(screen.getByText("f2.dat")).toBeInTheDocument();
});

test("zero-exposure sample renders disabled, no checkbox", () => {
  const empty: PickerSampleRow = {
    ...baseRow, indexing_exposure_id: null, all_exposures: [],
  };
  render(
    <SamplePickerRow
      row={empty} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  const row = screen.getByTestId("sample-picker-row");
  expect(row).toHaveAttribute("data-disabled", "true");
  expect(row).not.toHaveAttribute("data-exposure-id");
  expect(screen.queryByRole("checkbox")).toBeNull();
});

test("override radio fires onOverrideChange with selected exposure id", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={true} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={handle}
    />,
  );
  fireEvent.click(screen.getByTestId("sample-picker-row-caret"));
  fireEvent.click(screen.getByLabelText(/f2\.dat/));
  expect(handle).toHaveBeenCalledWith(101);
});

test("data-exposure-id reflects resolved id (default → indexing)", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={() => {}}
      overrideExposureId={undefined} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByTestId("sample-picker-row")).toHaveAttribute("data-exposure-id", "100");
});

test("data-exposure-id reflects override when set", () => {
  render(
    <SamplePickerRow
      row={baseRow} checked={true} onCheckedChange={() => {}}
      overrideExposureId={101} onOverrideChange={() => {}}
    />,
  );
  expect(screen.getByTestId("sample-picker-row")).toHaveAttribute("data-exposure-id", "101");
});

test("alreadyAdded renders locked row with hint", () => {
  const handle = vi.fn();
  render(
    <SamplePickerRow
      row={baseRow} checked={false} onCheckedChange={handle}
      overrideExposureId={undefined} onOverrideChange={() => {}}
      alreadyAdded={true}
    />,
  );
  const row = screen.getByTestId("sample-picker-row");
  expect(row).toHaveAttribute("data-locked", "true");
  expect(screen.getByText("already added")).toBeInTheDocument();
  const checkbox = screen.getByTestId("sample-picker-row-checkbox") as HTMLInputElement;
  expect(checkbox.disabled).toBe(true);
  expect(checkbox.checked).toBe(true);
  // Clicking the disabled checkbox must not fire onCheckedChange.
  fireEvent.click(checkbox);
  expect(handle).not.toHaveBeenCalled();
  // Caret still renders (user can inspect exposures).
  expect(screen.getByTestId("sample-picker-row-caret")).toBeInTheDocument();
});
