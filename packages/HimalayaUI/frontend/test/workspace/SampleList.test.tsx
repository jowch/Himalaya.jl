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
