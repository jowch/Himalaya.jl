import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { SampleFold } from "../src/print/components/SampleFold";
import type { LoadSample } from "../src/api";

function sample(over: Partial<LoadSample> = {}): LoadSample {
  return {
    sample_id: 10, name: "HA85 (S01P15)", slot_index: 15,
    grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
    exposures: [
      { id: 100, filename: "a1.tif", horizontal_position: 8, timestamp: "10:00" },
      { id: 101, filename: "a2.tif", horizontal_position: 36, timestamp: "10:01" },
    ],
    ...over,
  };
}

const noop = () => {};
const baseProps = {
  open: true, selected: false,
  onToggleOpen: noop, onToggleSelect: noop, onRename: noop,
  onSplit: noop, onMerge: noop, onDismissFlag: noop, onMoveExposure: noop,
  thumbSrcFor: () => null,
};

describe("SampleFold", () => {
  it("renders name + exposure count and the leaves when open", () => {
    render(<SampleFold sample={sample()} {...baseProps} />);
    expect(screen.getByText("HA85 (S01P15)")).toBeInTheDocument();
    expect(screen.getByText(/2 exposures/)).toBeInTheDocument();
    expect(screen.getAllByTestId("exposure-leaf")).toHaveLength(2);
  });

  it("checkbox toggles selection", () => {
    const onToggleSelect = vi.fn();
    render(<SampleFold sample={sample()} {...baseProps} onToggleSelect={onToggleSelect} />);
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onToggleSelect).toHaveBeenCalledWith(10);
  });

  it("merge flag shows the merge prompt with Merge / Keep separate", () => {
    const onMerge = vi.fn(); const onDismiss = vi.fn();
    render(<SampleFold sample={sample({ flag: { kind: "merge", merge_with_sample_id: 5, merge_with_label: "HA85 (S01P08)" } })}
      {...baseProps} onMerge={onMerge} onDismissFlag={onDismiss} />);
    expect(screen.getByText(/HA85 \(S01P08\)/)).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /merge/i }));
    expect(onMerge).toHaveBeenCalledWith(10, 5);
    fireEvent.click(screen.getByRole("button", { name: /keep separate/i }));
    expect(onDismiss).toHaveBeenCalledWith(10);
  });

  it("split flag shows a divider with the position jump before the split index", () => {
    render(<SampleFold sample={sample({ flag: { kind: "split", split_at_index: 1, jump_from: 8, jump_to: 36 } })} {...baseProps} />);
    const divider = screen.getByTestId("split-divider");
    expect(divider.textContent).toMatch(/8\.0 → 36\.0/);
  });

  describe("inline rename", () => {
    it("clicking Rename activates an inline Input with the current name", () => {
      render(<SampleFold sample={sample()} {...baseProps} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      expect(input).toBeInTheDocument();
      expect(input.value).toBe("HA85 (S01P15)");
    });

    it("Enter commits the trimmed new name and calls onRename", () => {
      const onRename = vi.fn();
      render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      fireEvent.change(input, { target: { value: "  New Name  " } });
      fireEvent.keyDown(input, { key: "Enter" });
      expect(onRename).toHaveBeenCalledWith(10, "New Name");
    });

    it("blur commits the trimmed new name and calls onRename", () => {
      const onRename = vi.fn();
      render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      fireEvent.change(input, { target: { value: "Blurred Name" } });
      fireEvent.blur(input);
      expect(onRename).toHaveBeenCalledWith(10, "Blurred Name");
    });

    it("Escape cancels rename without calling onRename", () => {
      const onRename = vi.fn();
      render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      fireEvent.change(input, { target: { value: "Should Not Save" } });
      fireEvent.keyDown(input, { key: "Escape" });
      expect(onRename).not.toHaveBeenCalled();
      // Back to static label
      expect(screen.getByText("HA85 (S01P15)")).toBeInTheDocument();
    });

    it("committing the same name (unchanged) does not call onRename", () => {
      const onRename = vi.fn();
      render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      // No change to value — commit via blur
      fireEvent.blur(input);
      expect(onRename).not.toHaveBeenCalled();
    });

    it("committing an empty name does not call onRename", () => {
      const onRename = vi.fn();
      render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
      fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
      const input = screen.getByTestId("sample-rename-input").querySelector("input")!;
      fireEvent.change(input, { target: { value: "   " } });
      fireEvent.blur(input);
      expect(onRename).not.toHaveBeenCalled();
    });
  });
});
