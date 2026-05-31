import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagEditor } from "../../src/print/ui/TagEditor";

describe("<TagEditor>", () => {
  it("commits a key-only tag (no value property) on Enter", () => {
    const onCommit = vi.fn();
    render(<TagEditor onCommit={onCommit} />);
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.change(keyInput, { target: { value: "lipid" } });
    fireEvent.keyDown(keyInput, { key: "Enter" });
    // exact object: a stray `value: undefined` would fail this assertion
    expect(onCommit).toHaveBeenCalledWith({ key: "lipid" });
  });

  it("commits a key+value tag on Enter", () => {
    const onCommit = vi.fn();
    render(<TagEditor onCommit={onCommit} />);
    fireEvent.change(screen.getByPlaceholderText("key"), {
      target: { value: "temperature" },
    });
    fireEvent.change(screen.getByPlaceholderText("value (optional)"), {
      target: { value: "37C" },
    });
    fireEvent.keyDown(screen.getByPlaceholderText("value (optional)"), {
      key: "Enter",
    });
    expect(onCommit).toHaveBeenCalledWith({ key: "temperature", value: "37C" });
  });

  it("does not commit when the key is empty", () => {
    const onCommit = vi.fn();
    render(<TagEditor onCommit={onCommit} />);
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.keyDown(keyInput, { key: "Enter" });
    expect(onCommit).not.toHaveBeenCalled();
  });

  it("calls onCancel on Escape", () => {
    const onCancel = vi.fn();
    render(<TagEditor onCommit={() => {}} onCancel={onCancel} />);
    fireEvent.keyDown(screen.getByPlaceholderText("key"), { key: "Escape" });
    expect(onCancel).toHaveBeenCalled();
  });
});
