import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SampleList } from "../../src/workspace/ui/SampleList";
import type { CorpusSample } from "../../src/api";

const sample = (id: number, name: string) => ({ id, name } as unknown as CorpusSample);

describe("SampleList filter keys", () => {
  it("leaves Enter to an input method that is still composing", () => {
    const onPick = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001")]} currentId={undefined} filter="JC" onFilter={() => {}} onPick={onPick} />,
    );
    const box = screen.getByLabelText("Filter samples");
    fireEvent.keyDown(box, { key: "Enter", isComposing: true });
    expect(onPick).not.toHaveBeenCalled();
    fireEvent.keyDown(box, { key: "Enter" });
    expect(onPick).toHaveBeenCalledWith(1);
  });
});

describe("SampleList weights", () => {
  it("sets only the current sample's name in semibold, so the cursor row stands out", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002")]} currentId={2} filter="" onFilter={() => {}} onPick={() => {}} />,
    );
    expect(screen.getByRole("button", { name: /JC002/ }).className).toMatch(/\bfont-semibold\b/);
    expect(screen.getByRole("button", { name: /JC001/ }).className).not.toMatch(/\bfont-semibold\b/);
  });
});

describe("SampleList comparisons", () => {
  it("marks compared rows", () => {
    Element.prototype.scrollIntoView = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001"), sample(2, "JC002")]} currentId={1} filter="" onFilter={() => {}} onPick={() => {}} compareColors={new Map([[2, "red"]])} />,
    );
    expect(screen.getByRole("button", { name: /JC002/ })).toHaveAttribute("data-compared");
    expect(screen.getByRole("button", { name: /JC001/ })).not.toHaveAttribute("data-compared");
  });

  it("a modifier click toggles the comparison instead of opening the sample", () => {
    const onPick = vi.fn();
    const onToggleCompare = vi.fn();
    render(
      <SampleList samples={[sample(1, "JC001")]} currentId={undefined} filter="" onFilter={() => {}} onPick={onPick} onToggleCompare={onToggleCompare} />,
    );
    const row = screen.getByRole("button", { name: /JC001/ });
    fireEvent.click(row, { metaKey: true });
    fireEvent.click(row, { ctrlKey: true });
    expect(onToggleCompare).toHaveBeenCalledTimes(2);
    expect(onPick).not.toHaveBeenCalled();
    fireEvent.click(row);
    expect(onPick).toHaveBeenCalledWith(1);
  });
});
