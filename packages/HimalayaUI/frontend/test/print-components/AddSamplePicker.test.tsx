import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { AddSamplePicker } from "../../src/print/components/AddSamplePicker";
import type { CorpusSample, Experiment } from "../../src/api";

function sample(id: number, name: string, experiment_id: number): CorpusSample {
  return {
    id,
    experiment_id,
    name,
    display_name: name,
    notes: null,
    tags: [],
    q_units: "A^-1",
  } as CorpusSample;
}
function exp(id: number, name: string): Experiment {
  return {
    id,
    name,
    path: `/x/${id}`,
    data_dir: "",
    analysis_dir: "",
    manifest_path: null,
    created_at: "2026-01-01",
    q_units: "A^-1",
  } as Experiment;
}

const EXPERIMENTS = [exp(1, "SSRL April"), exp(2, "APS March")];
const OPTIONS = [
  sample(10, "LL37 Only", 1),
  sample(11, "1-1 + LL37 1:1", 1),
  sample(20, "LL37 Only", 2), // duplicate display name, different experiment
];

function setup(over: Partial<React.ComponentProps<typeof AddSamplePicker>> = {}) {
  const onAdd = vi.fn();
  render(
    <AddSamplePicker
      options={OPTIONS}
      experiments={EXPERIMENTS}
      onAdd={onAdd}
      {...over}
    />,
  );
  return { onAdd };
}

function open(): void {
  fireEvent.click(screen.getByTestId("builder-add-sample"));
}

describe("<AddSamplePicker>", () => {
  it("opens on the trigger and focuses the search field (search-first)", () => {
    setup();
    expect(screen.queryByTestId("add-sample-listbox")).not.toBeInTheDocument();
    open();
    expect(screen.getByTestId("add-sample-listbox")).toBeInTheDocument();
    expect(screen.getByRole("combobox")).toHaveFocus();
  });

  it("groups options under experiment headers and carries each smp_{id}", () => {
    setup();
    open();
    expect(screen.getByRole("group", { name: "SSRL April" })).toBeInTheDocument();
    expect(screen.getByRole("group", { name: "APS March" })).toBeInTheDocument();
    // Every option row shows its mono id (disambiguates duplicate names).
    expect(within(screen.getByTestId("add-opt-10")).getByText("smp_10")).toBeInTheDocument();
    expect(within(screen.getByTestId("add-opt-20")).getByText("smp_20")).toBeInTheDocument();
  });

  it("filters by sample name", () => {
    setup();
    open();
    fireEvent.change(screen.getByRole("combobox"), { target: { value: "1:1" } });
    expect(screen.getByTestId("add-opt-11")).toBeInTheDocument();
    expect(screen.queryByTestId("add-opt-10")).not.toBeInTheDocument();
  });

  it("filters by smp_{id}", () => {
    setup();
    open();
    fireEvent.change(screen.getByRole("combobox"), { target: { value: "smp_20" } });
    expect(screen.getByTestId("add-opt-20")).toBeInTheDocument();
    expect(screen.queryByTestId("add-opt-10")).not.toBeInTheDocument();
  });

  it("filters by experiment name", () => {
    setup();
    open();
    fireEvent.change(screen.getByRole("combobox"), { target: { value: "APS" } });
    expect(screen.getByTestId("add-opt-20")).toBeInTheDocument();
    expect(screen.queryByTestId("add-opt-11")).not.toBeInTheDocument();
  });

  it("clicking an option adds it and keeps the picker open for multi-add", () => {
    const { onAdd } = setup();
    open();
    fireEvent.click(screen.getByTestId("add-opt-11"));
    expect(onAdd).toHaveBeenCalledWith(11);
    // Still open (the recipe-shrinks-options behaviour lives in the parent).
    expect(screen.getByTestId("add-sample-listbox")).toBeInTheDocument();
  });

  it("ArrowDown + Enter adds the active option from the keyboard (focus stays in the field)", () => {
    const { onAdd } = setup();
    open();
    const input = screen.getByRole("combobox");
    // Active defaults to the first option (smp_10); one ArrowDown → smp_11.
    fireEvent.keyDown(input, { key: "ArrowDown" });
    expect(input).toHaveAttribute("aria-activedescendant", "add-opt-11");
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onAdd).toHaveBeenCalledWith(11);
  });

  it("shows an honest empty state when nothing matches", () => {
    setup();
    open();
    fireEvent.change(screen.getByRole("combobox"), { target: { value: "zzz" } });
    expect(screen.getByText(/no samples match/i)).toBeInTheDocument();
    expect(screen.queryByTestId("add-opt-10")).not.toBeInTheDocument();
  });
});
