import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SampleList } from "../../src/ui/SampleList";
import type { CorpusSample } from "../../src/api";

const sample = (id: number, name: string) => ({ id, name } as unknown as CorpusSample);
const samples = [sample(1, "JC001"), sample(2, "JC002"), sample(3, "JC003")];
const list = (currentId: number | undefined) => (
  <SampleList samples={samples} currentId={currentId} filter="" onFilter={() => {}} onPick={() => {}} />
);
const row = (name: string) => screen.getByRole("button", { name: new RegExp(name) });

beforeEach(() => {
  Element.prototype.scrollIntoView = vi.fn();
});

describe("SampleList focus", () => {
  it("makes only the current row a Tab stop", () => {
    render(list(2));
    expect(row("JC001")).toHaveAttribute("tabindex", "-1");
    expect(row("JC002")).toHaveAttribute("tabindex", "0");
    expect(row("JC003")).toHaveAttribute("tabindex", "-1");
  });

  it("keeps the first row a Tab stop when no sample is current", () => {
    render(list(undefined));
    expect(row("JC001")).toHaveAttribute("tabindex", "0");
    expect(row("JC002")).toHaveAttribute("tabindex", "-1");
  });

  it("does not leave focus on a row picked with the mouse", () => {
    render(list(1));
    const r = row("JC002");
    r.focus();
    fireEvent.click(r, { detail: 1 });
    expect(r).not.toHaveFocus();
  });

  it("keeps focus on a row picked from the keyboard", () => {
    render(list(1));
    const r = row("JC002");
    r.focus();
    fireEvent.click(r, { detail: 0 });
    expect(r).toHaveFocus();
  });

  it("moves focus to the new current row when the selection moves while a row has focus", () => {
    const { rerender } = render(list(1));
    row("JC001").focus();
    rerender(list(2));
    expect(row("JC002")).toHaveFocus();
  });

  it("leaves focus alone when no row has it", () => {
    const { rerender } = render(<>{list(1)}<input aria-label="elsewhere" /></>);
    const other = screen.getByLabelText("elsewhere");
    other.focus();
    rerender(<>{list(2)}<input aria-label="elsewhere" /></>);
    expect(other).toHaveFocus();
  });
});
