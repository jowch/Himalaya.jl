import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { ComparePicker } from "../../src/workspace/ui/ComparePicker";
import { STANDARD_LABELS } from "../../src/lib/labels";
import type { CorpusSample } from "../../src/api";

const sample = (id: number, name: string) => ({ id, name, experiment_id: 1, notes: null, q_units: "A-1" }) as unknown as CorpusSample;

describe("ComparePicker", () => {
  it("keeps the highlighted sample when the list reorders under it", () => {
    const onPick = vi.fn();
    const a = sample(1, "A"), b = sample(2, "B"), c = sample(3, "C");
    const { rerender } = render(<ComparePicker samples={[a, b, c]} onPick={onPick} onClose={() => {}} />);
    const input = screen.getByLabelText("Find a sample to compare");
    fireEvent.keyDown(input, { key: "ArrowDown" });
    rerender(<ComparePicker samples={[c, a, b]} onPick={onPick} onClose={() => {}} />);
    expect(screen.getByRole("button", { current: true })).toHaveTextContent("B");
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onPick).toHaveBeenCalledWith(2);
  });

  it("reads a row's note as separate from its name", () => {
    render(<ComparePicker samples={[sample(1, "A3")]} notes={new Map([[1, "Buffer"]])} onPick={() => {}} onClose={() => {}} />);
    expect(within(screen.getByRole("dialog")).getByRole("button", { name: "A3, Buffer" })).toBeInTheDocument();
  });

  it("uses the same role names as the backend's ingest suggestions", () => {
    const jl = readFileSync(resolve(__dirname, "../../../src/label_suggestions.jl"), "utf8");
    const roles = [...jl.matchAll(/^\s+"([^"]+)"\s+=>/gm)].map((m) => m[1]);
    expect([...roles].sort()).toEqual([...STANDARD_LABELS].sort());
  });
});
