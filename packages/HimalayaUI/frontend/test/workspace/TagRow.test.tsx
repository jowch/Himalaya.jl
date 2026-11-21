import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import type { SampleLabel } from "../../src/api";
import { TagRow } from "../../src/workspace/ui/TagRow";

function tag(id: number, label: string, state: SampleLabel["state"] = "confirmed"): SampleLabel {
  return { id, sample_id: 1, label, state, source: "user" };
}

const NO_COLORS = { colors: [], coloring: undefined, onOpenColor: () => {}, onCloseColor: () => {}, onSetColor: () => {} };

const TAGS = [tag(1, "37C"), tag(2, "heated"), tag(3, "retake"), tag(4, "buffer", "suggested"), tag(5, "dilute")];

function renderRow(onAdd = vi.fn(), onRemove = vi.fn()): void {
  render(
    <TagRow
      tags={TAGS}
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

describe("TagRow overflow", () => {
  beforeEach(() => {
    vi.spyOn(HTMLElement.prototype, "offsetWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.parentElement?.dataset.testid === "workspace-tag-ruler" ? 60 : 0;
    });
    vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.dataset.testid === "workspace-tag-chips" ? 200 : 0;
    });
  });
  afterEach(() => vi.restoreAllMocks());

  it("shows the tags that fit and a +N chip that lists the rest, which stay editable", () => {
    const onAdd = vi.fn();
    const onRemove = vi.fn();
    renderRow(onAdd, onRemove);
    const chips = screen.getByTestId("workspace-tag-chips");
    expect([...chips.querySelectorAll("[data-tag]")].map((c) => c.getAttribute("data-tag"))).toEqual(["37C", "heated"]);
    fireEvent.click(screen.getByRole("button", { name: "3 more tags" }));
    const more = screen.getByRole("dialog", { name: "More tags" });
    expect([...more.querySelectorAll("[data-tag]")].map((c) => c.getAttribute("data-tag"))).toEqual(["retake", "buffer", "dilute"]);
    fireEvent.click(within(more).getByRole("button", { name: "Confirm tag Buffer" }));
    fireEvent.click(within(more).getByRole("button", { name: "Remove tag dilute" }));
    expect(onAdd).toHaveBeenCalledWith("buffer");
    expect(onRemove).toHaveBeenCalledWith("dilute");
    fireEvent.keyDown(more, { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "More tags" })).toBeNull();
  });

  it("reads a suggestion as suggested when it can't be confirmed", () => {
    render(
      <TagRow tags={[tag(4, "buffer", "suggested")]} vocabulary={[]} editable={false} adding={false}
        onOpenAdd={() => {}} onCloseAdd={() => {}} onAdd={() => {}} onRemove={() => {}} {...NO_COLORS} />,
    );
    expect(within(screen.getByTestId("workspace-tag-chips")).getByText(", suggested")).toBeInTheDocument();
  });
});

describe("TagRow colours", () => {
  it("right-clicking a tag in the +N list closes the list and opens that tag's colour menu", () => {
    vi.spyOn(HTMLElement.prototype, "offsetWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.parentElement?.dataset.testid === "workspace-tag-ruler" ? 60 : 0;
    });
    vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockImplementation(function (this: HTMLElement) {
      return this.dataset.testid === "workspace-tag-chips" ? 200 : 0;
    });
    const onOpenColor = vi.fn();
    render(
      <TagRow tags={TAGS} vocabulary={[]} editable adding={false} onOpenAdd={() => {}} onCloseAdd={() => {}}
        onAdd={() => {}} onRemove={() => {}} {...NO_COLORS} onOpenColor={onOpenColor} />,
    );
    fireEvent.click(screen.getByRole("button", { name: "3 more tags" }));
    fireEvent.contextMenu(screen.getByRole("dialog", { name: "More tags" }).querySelector("[data-tag='dilute']")!);
    expect(onOpenColor).toHaveBeenCalledWith("dilute");
    expect(screen.queryByRole("dialog", { name: "More tags" })).toBeNull();
    vi.restoreAllMocks();
  });
});
