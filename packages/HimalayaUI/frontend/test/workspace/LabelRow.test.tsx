import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import type { SampleLabel } from "../../src/api";
import { LabelRow } from "../../src/ui/LabelRow";

function chip(id: number, label: string, state: SampleLabel["state"] = "confirmed"): SampleLabel {
  return { id, sample_id: 1, label, state, source: "user" };
}

const NO_COLORS = { colors: [], coloring: undefined, onOpenColor: () => {}, onCloseColor: () => {}, onSetColor: () => {} };

const LABELS = [chip(1, "37C"), chip(2, "heated"), chip(3, "retake"), chip(4, "buffer", "suggested"), chip(5, "dilute")];

function renderRow(onAdd = vi.fn(), onRemove = vi.fn()): void {
  render(
    <LabelRow
      labels={LABELS}
      vocabulary={[]}
      editable
      adding={false}
      onOpenAdd={() => {}}
      onCloseAdd={() => {}}
      onAdd={onAdd}
      onRemove={onRemove}
      {...NO_COLORS}
    />,
  );
}

describe("LabelRow overflow", () => {
  beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.parentElement?.dataset.testid === "workspace-label-ruler" ? 60 : 0;
    });
    vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.dataset.testid === "workspace-label-chips" ? 200 : 0;
    });
  });
  afterEach(() => vi.restoreAllMocks());

  it("shows the labels that fit and a +N chip that lists the rest, which stay editable", () => {
    const onAdd = vi.fn();
    const onRemove = vi.fn();
    renderRow(onAdd, onRemove);
    const chips = screen.getByTestId("workspace-label-chips");
    expect([...chips.querySelectorAll("[data-label]")].map((c) => c.getAttribute("data-label"))).toEqual(["37C", "heated"]);
    fireEvent.click(screen.getByRole("button", { name: "3 more labels" }));
    const more = screen.getByRole("dialog", { name: "More labels" });
    expect([...more.querySelectorAll("[data-label]")].map((c) => c.getAttribute("data-label"))).toEqual(["retake", "buffer", "dilute"]);
    fireEvent.click(within(more).getByRole("button", { name: "Confirm label Buffer" }));
    fireEvent.click(within(more).getByRole("button", { name: "Remove label dilute" }));
    expect(onAdd).toHaveBeenCalledWith("buffer");
    expect(onRemove).toHaveBeenCalledWith("dilute");
    fireEvent.keyDown(more, { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "More labels" })).toBeNull();
  });

  it("reads a suggestion as suggested when it can't be confirmed", () => {
    render(
      <LabelRow labels={[chip(4, "buffer", "suggested")]} vocabulary={[]} editable={false} adding={false}
        onOpenAdd={() => {}} onCloseAdd={() => {}} onAdd={() => {}} onRemove={() => {}} {...NO_COLORS} />,
    );
    expect(within(screen.getByTestId("workspace-label-chips")).getByText(", suggested")).toBeInTheDocument();
  });
});

describe("LabelRow colours", () => {
  afterEach(() => vi.restoreAllMocks());

  function colorRow(over: Partial<Parameters<typeof LabelRow>[0]> = {}) {
    return (
      <LabelRow labels={LABELS} vocabulary={[]} editable adding={false} onOpenAdd={() => {}} onCloseAdd={() => {}}
        onAdd={() => {}} onRemove={() => {}} {...NO_COLORS} {...over} />
    );
  }

  it("right-clicking a label in the +N list closes the list and opens that label's colour menu", () => {
    vi.spyOn(HTMLElement.prototype, "offsetWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.parentElement?.dataset.testid === "workspace-label-ruler" ? 60 : 0;
    });
    vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.dataset.testid === "workspace-label-chips" ? 200 : 0;
    });
    const onOpenColor = vi.fn();
    render(colorRow({ onOpenColor }));
    fireEvent.click(screen.getByRole("button", { name: "3 more labels" }));
    fireEvent.contextMenu(screen.getByRole("dialog", { name: "More labels" }).querySelector("[data-label='dilute']")!);
    expect(onOpenColor).toHaveBeenCalledWith("dilute");
    expect(screen.queryByRole("dialog", { name: "More labels" })).toBeNull();
  });

  it("picking a label's default colour clears the override instead of storing it", () => {
    const onSetColor = vi.fn();
    render(colorRow({ coloring: "heated", colors: [{ label: "heated", color: "teal" }], onSetColor }));
    fireEvent.click(screen.getByRole("radio", { name: "Violet" }));
    fireEvent.click(screen.getByRole("radio", { name: "Rose" }));
    expect(onSetColor.mock.calls).toEqual([["heated", null], ["heated", "rose"]]);
  });

  it("closing hands focus back to what opened the menu, never to a label's remove button", () => {
    const { rerender } = render(<><button type="button">opener</button>{colorRow()}</>);
    screen.getByRole("button", { name: "opener" }).focus();
    rerender(<><button type="button">opener</button>{colorRow({ coloring: "37C" })}</>);
    fireEvent.keyDown(screen.getByRole("radio", { name: "Violet" }), { key: "Escape" });
    rerender(<><button type="button">opener</button>{colorRow({ coloring: undefined })}</>);
    expect(screen.getByRole("button", { name: "opener" })).toHaveFocus();

    screen.getByRole("button", { name: "Remove label 37C" }).focus();
    rerender(<><button type="button">opener</button>{colorRow({ coloring: "37C" })}</>);
    fireEvent.click(screen.getByRole("radio", { name: "Rose" }));
    rerender(<><button type="button">opener</button>{colorRow({ coloring: undefined })}</>);
    expect(screen.getByRole("button", { name: "Remove label 37C" })).not.toHaveFocus();
    expect(screen.getByRole("button", { name: "+ Label L" })).toHaveFocus();
  });

  it("falls back to + Label when the opener is gone, as after the palette closes", () => {
    const { rerender } = render(colorRow({ coloring: "37C" }));
    fireEvent.keyDown(screen.getByRole("radio", { name: "Violet" }), { key: "Escape" });
    rerender(colorRow({ coloring: undefined }));
    expect(screen.getByRole("button", { name: "+ Label L" })).toHaveFocus();
  });
});
