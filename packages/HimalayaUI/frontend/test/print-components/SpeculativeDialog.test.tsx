import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SpeculativeDialog } from "../../src/print/components/SpeculativeDialog";
import type { SpeculativeSnap } from "../../src/api";

const PEAKS = [
  { id: 1, q: 0.0301, source: "auto" },
  { id: 2, q: 0.0426, source: "auto" },
  { id: 3, q: 0.0521, source: "curation" },
];

// One anchor row, one row with a suggested peak, one with no peak.
const SNAP: SpeculativeSnap[] = [
  {
    ratio_position: 1,
    predicted_q: 0.0301,
    suggested_peak_id: 1,
    suggested_q: 0.0301,
    suggested_residual: 0,
    is_anchor: true,
  },
  {
    ratio_position: 2,
    predicted_q: 0.0426,
    suggested_peak_id: 2,
    suggested_q: 0.0428,
    suggested_residual: 0.0002,
    is_anchor: false,
  },
  {
    ratio_position: 3,
    predicted_q: 0.0521,
    suggested_peak_id: null,
    suggested_q: null,
    suggested_residual: null,
    is_anchor: false,
  },
];

function baseProps(over: Partial<React.ComponentProps<typeof SpeculativeDialog>> = {}) {
  return {
    open: true,
    onClose: vi.fn(),
    phases: ["Lamellar", "Hexagonal", "Pn3m"] as readonly string[],
    phase: "Lamellar",
    onPhaseChange: vi.fn(),
    peaks: PEAKS,
    anchorPeakId: 1 as number | undefined,
    onAnchorChange: vi.fn(),
    anchorRatio: 1,
    onAnchorRatioChange: vi.fn(),
    snap: SNAP,
    included: {} as Record<number, boolean>,
    onToggleIncluded: vi.fn(),
    onCreate: vi.fn(),
    ...over,
  };
}

describe("SpeculativeDialog", () => {
  it("renders nothing when open=false", () => {
    const { queryByText } = render(<SpeculativeDialog {...baseProps({ open: false })} />);
    expect(queryByText("Speculative index")).toBeNull();
  });

  it("renders title, phase/anchor selects, ratio input and snap rows when open", () => {
    const { getByText, getByTestId } = render(<SpeculativeDialog {...baseProps()} />);
    expect(getByText("Speculative index")).toBeTruthy();
    expect(getByTestId("spec-phase-select")).toBeTruthy();
    expect(getByTestId("spec-anchor-select")).toBeTruthy();
    expect(getByTestId("spec-snap-row-1")).toBeTruthy();
    expect(getByTestId("spec-snap-row-2")).toBeTruthy();
    expect(getByTestId("spec-snap-row-3")).toBeTruthy();
    // ratio input present (number)
    const inputs = document.querySelectorAll('input[type="number"]');
    expect(inputs.length).toBe(1);
  });

  it("fires onPhaseChange / onAnchorChange(number) / onAnchorRatioChange on edits", () => {
    const props = baseProps();
    const { getByTestId } = render(<SpeculativeDialog {...props} />);
    fireEvent.change(getByTestId("spec-phase-select"), { target: { value: "Pn3m" } });
    expect(props.onPhaseChange).toHaveBeenCalledWith("Pn3m");

    fireEvent.change(getByTestId("spec-anchor-select"), { target: { value: "2" } });
    expect(props.onAnchorChange).toHaveBeenCalledWith(2);

    const input = document.querySelector('input[type="number"]') as HTMLInputElement;
    fireEvent.change(input, { target: { value: "3" } });
    expect(props.onAnchorRatioChange).toHaveBeenCalledWith(3);
  });

  it("toggles a non-anchor has-peak row; disables anchor and no-peak checkboxes", () => {
    const props = baseProps();
    const { getByTestId } = render(<SpeculativeDialog {...props} />);
    const cb = (row: number) =>
      getByTestId(`spec-snap-row-${row}`).querySelector('input[type="checkbox"]') as HTMLInputElement;

    expect(cb(1).disabled).toBe(true); // anchor
    expect(cb(3).disabled).toBe(true); // no peak
    expect(cb(2).disabled).toBe(false);

    fireEvent.click(cb(2));
    expect(props.onToggleIncluded).toHaveBeenCalledWith(2);
  });

  it("disables Save when no anchor; enables and fires onCreate with an anchor", () => {
    const noAnchor = baseProps({ anchorPeakId: undefined });
    const { getByTestId, rerender } = render(<SpeculativeDialog {...noAnchor} />);
    expect((getByTestId("spec-save-button") as HTMLButtonElement).disabled).toBe(true);

    const withAnchor = baseProps();
    rerender(<SpeculativeDialog {...withAnchor} />);
    const save = getByTestId("spec-save-button") as HTMLButtonElement;
    expect(save.disabled).toBe(false);
    fireEvent.click(save);
    expect(withAnchor.onCreate).toHaveBeenCalledTimes(1);
  });

  it("shows Saving… and disables Save while saving", () => {
    const { getByTestId } = render(<SpeculativeDialog {...baseProps({ saving: true })} />);
    const save = getByTestId("spec-save-button") as HTMLButtonElement;
    expect(save.disabled).toBe(true);
    expect(save.textContent).toContain("Saving…");
  });

  it("renders the create error in a role=alert", () => {
    const { getByRole } = render(<SpeculativeDialog {...baseProps({ error: "boom" })} />);
    expect(getByRole("alert").textContent).toContain("boom");
  });
});
