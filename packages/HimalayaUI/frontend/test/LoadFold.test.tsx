import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { LoadFold } from "../src/print/components/LoadFold";
import type { Load } from "../src/api";

function load(over: Partial<Load> = {}): Load {
  return {
    load_id: 1, load_index: 1, session_id: null, start_time: "10:02", end_time: "10:38", frame_count: 96, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "merge", merge_with_sample_id: 5, merge_with_label: "C" }, exposures: [] },
    ],
    ...over,
  };
}

const cb = {
  onToggleLoad: () => {}, openSamples: new Set<number>(), selected: new Set<number>(),
  onToggleSampleOpen: () => {}, onToggleSelect: () => {}, onRename: () => {},
  onSplit: () => {}, onMerge: () => {}, onDismissFlag: () => {}, onMoveExposure: () => {},
  thumbSrcFor: () => null,
};

describe("LoadFold", () => {
  it("shows title and an attn status with the flagged count when not clean", () => {
    render(<LoadFold load={load()} open visibleSamples={load().samples} {...cb} />);
    expect(screen.getByText("Load 1")).toBeInTheDocument();
    expect(screen.getByText(/1 to check/)).toBeInTheDocument();
  });
  it("shows a clean status when no sample is flagged", () => {
    const clean = load({ samples: [load().samples[0]!] });
    render(<LoadFold load={clean} open visibleSamples={clean.samples} {...cb} />);
    expect(screen.getByText(/grouped cleanly/i)).toBeInTheDocument();
  });
  it("renders the subset hint when fewer samples are visible than total", () => {
    const l = load();
    render(<LoadFold load={l} open visibleSamples={[l.samples[1]!]} {...cb} />);
    expect(screen.getByText(/1 of 2 shown/)).toBeInTheDocument();
  });
  it("header click toggles the load", () => {
    const onToggleLoad = vi.fn();
    render(<LoadFold load={load()} open={false} visibleSamples={[]} {...cb} onToggleLoad={onToggleLoad} />);
    fireEvent.click(screen.getByRole("button", { name: /load 1/i }));
    expect(onToggleLoad).toHaveBeenCalledWith(1);
  });
});
