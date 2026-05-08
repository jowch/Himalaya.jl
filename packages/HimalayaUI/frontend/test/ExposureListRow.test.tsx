/**
 * ExposureListRow (Plan §Phase 5, Task 5.1).
 *
 * The row component shared by the picker modal and (in future) any list-mode
 * exposure surface. Renders exposure name + sample name + (truncated) sample
 * notes; takes a slot for trailing controls (checkbox / click button).
 */
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { ExposureListRow } from "../src/components/ExposureListRow";
import type { Exposure, Sample } from "../src/api";

function makeExposure(overrides: Partial<Exposure> = {}): Exposure {
  return {
    id: 100,
    sample_id: 10,
    filename: "JC001-120.dat",
    kind: "file",
    selected: false,
    status: "accepted",
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
    ...overrides,
  };
}

function makeSample(overrides: Partial<Sample> = {}): Sample {
  return {
    id: 10,
    experiment_id: 1,
    label: "JC001",
    name: "Sample A1",
    notes: "DOPC + 50% chol, 10 mM CaCl2, hydrated 24h",
    tags: [],
    ...overrides,
  };
}

describe("<ExposureListRow>", () => {
  it("renders exposure name (filename without .dat) and sample name", () => {
    render(
      <ExposureListRow
        exposure={makeExposure({ filename: "JC001-120.dat" })}
        sample={makeSample({ name: "Sample A1" })}
      />,
    );
    expect(screen.getByText("JC001-120")).toBeInTheDocument();
    expect(screen.getByText("Sample A1")).toBeInTheDocument();
  });

  it("falls back to sample.label when sample.name is missing", () => {
    render(
      <ExposureListRow
        exposure={makeExposure()}
        sample={makeSample({ name: null, label: "JC042" })}
      />,
    );
    expect(screen.getByText("JC042")).toBeInTheDocument();
  });

  it("renders notes (truncated) and exposes full text via title attribute", () => {
    const longNotes =
      "this is a very long sample note that should be truncated to about two lines visually but the full text remains available via the title attribute for hover";
    render(
      <ExposureListRow
        exposure={makeExposure()}
        sample={makeSample({ notes: longNotes })}
      />,
    );
    const notesEl = screen.getByTestId("exposure-list-row-notes");
    expect(notesEl).toHaveAttribute("title", longNotes);
    expect(notesEl).toHaveTextContent(longNotes);
  });

  it("renders no notes element when sample.notes is null", () => {
    render(
      <ExposureListRow
        exposure={makeExposure()}
        sample={makeSample({ notes: null })}
      />,
    );
    expect(screen.queryByTestId("exposure-list-row-notes")).toBeNull();
  });

  it("invokes onClick when the row is clicked (default action mode)", async () => {
    const user = userEvent.setup();
    const onClick = vi.fn();
    render(
      <ExposureListRow
        exposure={makeExposure({ id: 7 })}
        sample={makeSample()}
        onClick={onClick}
      />,
    );
    await user.click(screen.getByTestId("exposure-list-row"));
    expect(onClick).toHaveBeenCalledTimes(1);
    expect(onClick).toHaveBeenCalledWith(7);
  });

  it("renders an actionElement slot at the trailing edge", () => {
    render(
      <ExposureListRow
        exposure={makeExposure()}
        sample={makeSample()}
        actionElement={<span data-testid="custom-action">X</span>}
      />,
    );
    expect(screen.getByTestId("custom-action")).toBeInTheDocument();
  });

  it("checkbox mode: renders a checkbox and toggles via onCheckedChange", async () => {
    const user = userEvent.setup();
    const onCheckedChange = vi.fn();
    render(
      <ExposureListRow
        exposure={makeExposure({ id: 9 })}
        sample={makeSample()}
        checked={false}
        onCheckedChange={onCheckedChange}
      />,
    );
    const checkbox = screen.getByTestId("exposure-list-row-checkbox");
    await user.click(checkbox);
    expect(onCheckedChange).toHaveBeenCalledTimes(1);
    expect(onCheckedChange).toHaveBeenCalledWith(true);
  });

  it("checkbox mode: when locked, renders disabled checkbox and surfaces data-locked", () => {
    const onCheckedChange = vi.fn();
    render(
      <ExposureListRow
        exposure={makeExposure({ id: 9 })}
        sample={makeSample()}
        checked
        onCheckedChange={onCheckedChange}
        locked
        lockedReason="already added"
      />,
    );
    const row = screen.getByTestId("exposure-list-row");
    expect(row).toHaveAttribute("data-locked", "true");
    const checkbox = screen.getByTestId(
      "exposure-list-row-checkbox",
    ) as HTMLInputElement;
    expect(checkbox.disabled).toBe(true);
    // Still shows a hint of why it's locked
    expect(screen.getByText(/already added/i)).toBeInTheDocument();
  });

  it("data attributes are stable: data-exposure-id and data-testid='exposure-list-row'", () => {
    render(
      <ExposureListRow
        exposure={makeExposure({ id: 42 })}
        sample={makeSample()}
      />,
    );
    const row = screen.getByTestId("exposure-list-row");
    expect(row.getAttribute("data-exposure-id")).toBe("42");
  });

  it("clicking a locked row does not fire onClick", async () => {
    const user = userEvent.setup();
    const onClick = vi.fn();
    render(
      <ExposureListRow
        exposure={makeExposure({ id: 7 })}
        sample={makeSample()}
        onClick={onClick}
        locked
      />,
    );
    await user.click(screen.getByTestId("exposure-list-row"));
    expect(onClick).not.toHaveBeenCalled();
  });

  it("renders radio when control='radio'", () => {
    render(
      <ExposureListRow
        exposure={makeExposure()} sample={makeSample()}
        control="radio"
        checked={false} onCheckedChange={() => {}}
      />,
    );
    const input = screen.getByTestId("exposure-list-row-checkbox") as HTMLInputElement;
    expect(input.type).toBe("radio");
  });

  it("renders checkbox when control unset (default)", () => {
    render(
      <ExposureListRow
        exposure={makeExposure()} sample={makeSample()}
        checked={false} onCheckedChange={() => {}}
      />,
    );
    const input = screen.getByTestId("exposure-list-row-checkbox") as HTMLInputElement;
    expect(input.type).toBe("checkbox");
  });
});
