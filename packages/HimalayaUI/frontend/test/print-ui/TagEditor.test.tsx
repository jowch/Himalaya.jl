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

  it("rejects an empty-key Add honestly: aria-invalid on the key field + inline error text", () => {
    const onCommit = vi.fn();
    render(<TagEditor onCommit={onCommit} />);
    // Quiet until the user actually tries to commit.
    expect(screen.queryByTestId("tag-editor-error")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(onCommit).not.toHaveBeenCalled();
    const keyInput = screen.getByPlaceholderText("key");
    expect(keyInput).toHaveAttribute("aria-invalid", "true");
    const error = screen.getByTestId("tag-editor-error");
    expect(error).toHaveTextContent("Enter a key first.");
    // The error is programmatically tied to the field, not just adjacent.
    expect(keyInput).toHaveAttribute("aria-describedby", error.getAttribute("id"));
  });

  it("Enter with an empty key raises the same inline error", () => {
    render(<TagEditor onCommit={() => {}} />);
    fireEvent.keyDown(screen.getByPlaceholderText("key"), { key: "Enter" });
    expect(screen.getByTestId("tag-editor-error")).toBeInTheDocument();
  });

  it("typing a key clears the error", () => {
    render(<TagEditor onCommit={() => {}} />);
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(screen.getByTestId("tag-editor-error")).toBeInTheDocument();
    fireEvent.change(keyInput, { target: { value: "l" } });
    expect(screen.queryByTestId("tag-editor-error")).toBeNull();
    expect(keyInput).not.toHaveAttribute("aria-invalid");
  });

  it("calls onCancel on Escape", () => {
    const onCancel = vi.fn();
    render(<TagEditor onCommit={() => {}} onCancel={onCancel} />);
    fireEvent.keyDown(screen.getByPlaceholderText("key"), { key: "Escape" });
    expect(onCancel).toHaveBeenCalled();
  });

  // ── existingKeys: single-valued-key rule ────────────────────────────────────

  it("rejects a key already in existingKeys: aria-invalid + inline error text, no commit", () => {
    const onCommit = vi.fn();
    render(
      <TagEditor
        onCommit={onCommit}
        existingKeys={["lipid", "temperature"]}
      />,
    );
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.change(keyInput, { target: { value: "lipid" } });
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(onCommit).not.toHaveBeenCalled();
    expect(keyInput).toHaveAttribute("aria-invalid", "true");
    const error = screen.getByTestId("tag-editor-error");
    expect(error).toHaveTextContent(/already has a tag with that key/);
    expect(keyInput).toHaveAttribute("aria-describedby", error.getAttribute("id"));
  });

  it("does NOT reject a key that is NOT in existingKeys", () => {
    const onCommit = vi.fn();
    render(
      <TagEditor onCommit={onCommit} existingKeys={["lipid"]} />,
    );
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.change(keyInput, { target: { value: "dose" } });
    fireEvent.keyDown(keyInput, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith({ key: "dose" });
  });

  it("typing a key clears the existingKeys duplicate error", () => {
    render(
      <TagEditor onCommit={() => {}} existingKeys={["lipid"]} />,
    );
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.change(keyInput, { target: { value: "lipid" } });
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(screen.getByTestId("tag-editor-error")).toBeInTheDocument();
    fireEvent.change(keyInput, { target: { value: "l" } });
    expect(screen.queryByTestId("tag-editor-error")).toBeNull();
    expect(keyInput).not.toHaveAttribute("aria-invalid");
  });

  it("without existingKeys, a key matching another tag commits normally", () => {
    // existingKeys is optional — omitting it skips the check (backward-compat)
    const onCommit = vi.fn();
    render(<TagEditor onCommit={onCommit} />);
    const keyInput = screen.getByPlaceholderText("key");
    fireEvent.change(keyInput, { target: { value: "lipid" } });
    fireEvent.keyDown(keyInput, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith({ key: "lipid" });
  });
});
