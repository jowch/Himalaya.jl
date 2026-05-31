import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagList } from "../../src/print/ui/TagList";

describe("TagList", () => {
  it("renders all tags as pills", () => {
    render(<TagList tags={["lipid-A", "37C", "run-3"]} />);
    expect(screen.getAllByTestId("tag-pill")).toHaveLength(3);
  });

  it("renders NO add button when onAdd is omitted", () => {
    render(<TagList tags={["lipid-A"]} />);
    expect(screen.queryByTestId("tag-add")).toBeNull();
  });

  it("renders an add button with aria-label when onAdd is provided and fires it on click", () => {
    const onAdd = vi.fn();
    render(<TagList tags={["lipid-A"]} onAdd={onAdd} />);
    const add = screen.getByTestId("tag-add");
    expect(add).toBeTruthy();
    expect(add.getAttribute("aria-label")).toBe("Add tag");
    fireEvent.click(add);
    expect(onAdd).toHaveBeenCalledTimes(1);
  });

  it("renders a remove button per pill when editable and onRemove are provided", () => {
    const onRemove = vi.fn();
    render(
      <TagList tags={["lipid-A", "37C"]} editable onRemove={onRemove} />,
    );
    expect(
      screen.getAllByRole("button", { name: "Remove tag" }),
    ).toHaveLength(2);
  });

  it("renders NO remove buttons when not editable", () => {
    const onRemove = vi.fn();
    render(<TagList tags={["lipid-A", "37C"]} onRemove={onRemove} />);
    expect(screen.queryAllByRole("button", { name: "Remove tag" })).toHaveLength(
      0,
    );
  });

  it("shows the add invite even when the tag list is empty", () => {
    const onAdd = vi.fn();
    render(<TagList tags={[]} onAdd={onAdd} />);
    expect(screen.getByTestId("tag-add")).toBeTruthy();
  });
});
