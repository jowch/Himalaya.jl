import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { LabelColorMenu } from "../../src/workspace/ui/LabelColorMenu";

function menu(onPick = vi.fn(), onClose = vi.fn()) {
  return <LabelColorMenu name="37C" current="violet" isDefault left={0} onPick={onPick} onClose={onClose} />;
}

describe("LabelColorMenu", () => {
  it("focuses the current colour once and keeps focus where the user moved it across rerenders", () => {
    const { rerender } = render(menu());
    expect(screen.getByRole("radio", { name: "Violet" })).toHaveFocus();
    fireEvent.keyDown(screen.getByRole("radio", { name: "Violet" }), { key: "ArrowRight" });
    expect(screen.getByRole("radio", { name: "Blue" })).toHaveFocus();
    rerender(menu(vi.fn(), vi.fn()));
    expect(screen.getByRole("radio", { name: "Blue" })).toHaveFocus();
  });

  it("is one tab stop whose arrows, Home and End move through the colours, wrapping", () => {
    render(menu());
    expect(screen.getAllByRole("radio").map((r) => r.tabIndex)).toEqual([0, -1, -1, -1, -1]);
    const group = screen.getByRole("radiogroup");
    fireEvent.keyDown(group, { key: "ArrowLeft" });
    expect(screen.getByRole("radio", { name: "Stone" })).toHaveFocus();
    fireEvent.keyDown(group, { key: "Home" });
    expect(screen.getByRole("radio", { name: "Violet" })).toHaveFocus();
    fireEvent.keyDown(group, { key: "End" });
    expect(screen.getByRole("radio", { name: "Stone" })).toHaveFocus();
  });

  it("Escape and a pick close it handing focus back; an outside click closes it without", () => {
    const onPick = vi.fn();
    const onClose = vi.fn();
    render(menu(onPick, onClose));
    fireEvent.keyDown(screen.getByRole("radio", { name: "Violet" }), { key: "Escape" });
    fireEvent.click(screen.getByRole("radio", { name: "Rose" }));
    fireEvent.mouseDown(document.body);
    expect(onPick).toHaveBeenCalledWith("rose");
    expect(onClose.mock.calls).toEqual([[true], [true], [false]]);
  });
});
