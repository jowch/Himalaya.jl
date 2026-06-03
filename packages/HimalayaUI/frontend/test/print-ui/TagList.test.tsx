import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagList } from "../../src/print/ui/TagList";
import type { Tag } from "../../src/print/ui/tag";

describe("TagList", () => {
  it("renders one TagPill per tag, including key+value tags", () => {
    const tags: Tag[] = [
      { key: "LL37" },
      { key: "temperature", value: "37C" },
      { key: "buffer", value: "PBS" },
    ];
    render(<TagList tags={tags} />);
    expect(screen.getAllByTestId("tag-pill")).toHaveLength(3);
    expect(
      screen.getAllByTestId("tag-pill")[0]!.getAttribute("data-size"),
    ).toBe("sm");
    expect(
      screen
        .getAllByTestId("tag-pill")
        .some(
          (p) =>
            p.querySelector('[data-role="tag-value"]')?.textContent === "37C",
        ),
    ).toBe(true);
  });

  it("renders NO add chip when onAdd is omitted", () => {
    render(<TagList tags={[{ key: "LL37" }]} />);
    expect(screen.queryByRole("button", { name: "Add" })).toBeNull();
  });

  it("shows the add chip when onAdd is provided; clicking it reveals the editor, and committing a key fires onAdd and hides the editor", () => {
    const onAdd = vi.fn();
    render(<TagList tags={[{ key: "LL37" }]} onAdd={onAdd} />);

    const add = screen.getByRole("button", { name: "Add" });
    expect(add.getAttribute("data-variant")).toBe("add");
    expect(screen.queryByTestId("tag-editor")).toBeNull();

    fireEvent.click(add);
    const editor = screen.getByTestId("tag-editor");
    expect(editor).toBeTruthy();

    const keyInput = editor.querySelector(
      'input[placeholder="key"]',
    ) as HTMLInputElement;
    fireEvent.change(keyInput, { target: { value: "pH" } });
    fireEvent.keyDown(keyInput, { key: "Enter" });

    expect(onAdd).toHaveBeenCalledTimes(1);
    expect(onAdd).toHaveBeenCalledWith({ key: "pH" });
    expect(screen.queryByTestId("tag-editor")).toBeNull();
  });

  it("renders a remove button per pill when editable and onRemove are provided, firing onRemove with the tag", () => {
    const onRemove = vi.fn();
    const tags: Tag[] = [{ key: "LL37" }, { key: "temperature", value: "37C" }];
    render(<TagList tags={tags} editable onRemove={onRemove} />);
    const buttons = screen.getAllByRole("button", { name: "Remove" });
    expect(buttons).toHaveLength(2);
    fireEvent.click(buttons[1]!);
    expect(onRemove).toHaveBeenCalledWith({ key: "temperature", value: "37C" });
  });

  it("renders NO remove buttons when not editable", () => {
    const onRemove = vi.fn();
    render(
      <TagList tags={[{ key: "LL37" }, { key: "x", value: "y" }]} onRemove={onRemove} />,
    );
    expect(screen.queryAllByRole("button", { name: "Remove" })).toHaveLength(0);
  });

  it("shows the add chip even when the tag list is empty (the invite)", () => {
    const onAdd = vi.fn();
    render(<TagList tags={[]} onAdd={onAdd} />);
    expect(screen.getByRole("button", { name: "Add" })).toBeTruthy();
  });

  it("does NOT flag persistent-add by default (hover-gated when tags present)", () => {
    render(<TagList tags={[{ key: "LL37" }]} onAdd={vi.fn()} />);
    expect(screen.getByTestId("tag-list")).not.toHaveAttribute(
      "data-persistent-add",
    );
  });

  it("persistentAdd keeps the add invite always shown alongside existing tags (the loupe rule)", () => {
    render(<TagList tags={[{ key: "LL37" }]} onAdd={vi.fn()} persistentAdd />);
    expect(screen.getByTestId("tag-list")).toHaveAttribute(
      "data-persistent-add",
      "true",
    );
    // The add affordance is still the functional dashed add chip.
    expect(screen.getByRole("button", { name: "Add" }).getAttribute("data-variant")).toBe(
      "add",
    );
  });

  describe("maxVisible overflow cap", () => {
    const FIVE: Tag[] = [
      { key: "LL37" },
      { key: "temperature", value: "37C" },
      { key: "buffer", value: "PBS" },
      { key: "lipid", value: "DOPC" },
      { key: "pH", value: "7.4" },
    ];

    it("caps at maxVisible pills plus one muted '+N' overflow trigger (NOT a tag-pill)", () => {
      render(<TagList tags={FIVE} maxVisible={2} />);
      expect(screen.getAllByTestId("tag-pill")).toHaveLength(2);
      const overflow = screen.getByTestId("tag-overflow");
      expect(overflow.textContent).toBe("+3 more");
      // The trigger is a real <button>, and it is NOT styled as a tag pill.
      expect(overflow.tagName).toBe("BUTTON");
      expect(overflow.getAttribute("data-testid")).not.toBe("tag-pill");
    });

    it("clicking the overflow trigger opens a popover with the hidden tags as real pills", () => {
      render(<TagList tags={FIVE} maxVisible={2} />);
      // Closed: only the 2 visible pills, no popover.
      expect(screen.getAllByTestId("tag-pill")).toHaveLength(2);
      expect(screen.queryByTestId("popover")).toBeNull();

      fireEvent.click(screen.getByTestId("tag-overflow"));
      const panel = screen.getByTestId("popover");
      expect(panel).toBeTruthy();
      // The 3 hidden tags now render as real TagPills inside the popover.
      const pillsInPanel = panel.querySelectorAll('[data-testid="tag-pill"]');
      expect(pillsInPanel).toHaveLength(3);
      expect(panel.textContent).toContain("buffer");
      expect(panel.textContent).toContain("PBS");
      expect(panel.textContent).toContain("7.4");
    });

    it("wires onRemove on the hidden popover pills when editable", () => {
      const onRemove = vi.fn();
      render(<TagList tags={FIVE} maxVisible={2} editable onRemove={onRemove} />);
      fireEvent.click(screen.getByTestId("tag-overflow"));
      const panel = screen.getByTestId("popover");
      const removes = panel.querySelectorAll('button[aria-label="Remove"]');
      expect(removes).toHaveLength(3);
      fireEvent.click(removes[0]!);
      expect(onRemove).toHaveBeenCalledWith({ key: "buffer", value: "PBS" });
    });

    it("renders all pills and no overflow chip when maxVisible is unset", () => {
      render(<TagList tags={FIVE} />);
      expect(screen.getAllByTestId("tag-pill")).toHaveLength(5);
      expect(screen.queryByTestId("tag-overflow")).toBeNull();
    });

    it("renders no overflow chip when tags.length <= maxVisible", () => {
      render(<TagList tags={FIVE.slice(0, 2)} maxVisible={2} />);
      expect(screen.getAllByTestId("tag-pill")).toHaveLength(2);
      expect(screen.queryByTestId("tag-overflow")).toBeNull();
    });
  });
});
