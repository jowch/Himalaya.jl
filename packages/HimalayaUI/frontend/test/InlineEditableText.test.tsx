import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { InlineEditableText } from "../src/components/InlineEditableText";

describe("InlineEditableText — Compare UX C-3", () => {
  it("renders text by default, switches to input on click", () => {
    render(<InlineEditableText value="hello" onCommit={() => {}} />);
    const text = screen.getByText("hello");
    fireEvent.click(text);
    expect(screen.getByRole("textbox")).toHaveValue("hello");
  });

  it("commits on Enter", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    const input = screen.getByRole("textbox");
    fireEvent.change(input, { target: { value: "world" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("world");
  });

  it("cancels on Esc", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "world" } });
    fireEvent.keyDown(screen.getByRole("textbox"), { key: "Escape" });
    expect(onCommit).not.toHaveBeenCalled();
    expect(screen.getByText("hello")).toBeInTheDocument();
  });

  it("commits on blur", () => {
    const onCommit = vi.fn();
    render(<InlineEditableText value="hello" onCommit={onCommit} />);
    fireEvent.click(screen.getByText("hello"));
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "world" } });
    fireEvent.blur(screen.getByRole("textbox"));
    expect(onCommit).toHaveBeenCalledWith("world");
  });

  it("carries data-interactable='edit'", () => {
    render(<InlineEditableText value="hello" onCommit={() => {}} />);
    expect(screen.getByText("hello").closest("[data-interactable]"))
      .toHaveAttribute("data-interactable", "edit");
  });

  it("does not become editable when readOnly", () => {
    render(<InlineEditableText value="locked" onCommit={vi.fn()} readOnly />);
    fireEvent.click(screen.getByText("locked"));
    expect(screen.queryByRole("textbox")).toBeNull();
    expect(screen.getByText("locked").closest("[data-interactable]")).toBeNull();
  });
});
