import { render, screen, fireEvent, within, act } from "@testing-library/react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { createRef } from "react";
import { ManageTagsModal } from "../../src/print/components/ManageTagsModal";
import type { LoupeTag } from "../../src/print/pages/loupeAdapters";

const KEY_OPTIONS = [
  { label: "lipid", count: 4 },
  { label: "temperature", count: 2 },
  { label: "dose", count: 3 },
];

function valueOptions(key: string) {
  const map: Record<string, { label: string; count?: number }[]> = {
    lipid: [{ label: "LL37", count: 3 }, { label: "DMPC", count: 1 }],
    temperature: [{ label: "37C", count: 2 }, { label: "25C", count: 1 }],
    dose: [{ label: "10", count: 2 }],
  };
  return map[key] ?? [];
}

const TAGS: LoupeTag[] = [
  { id: 1, key: "lipid", value: "LL37", source: "manual" },
  { id: 2, key: "temperature", value: "37C", source: "manual" },
];

function baseProps(over: Partial<Parameters<typeof ManageTagsModal>[0]> = {}) {
  return {
    open: true,
    sampleName: "JC042 — LL37",
    tags: TAGS,
    keyOptions: KEY_OPTIONS,
    valueOptionsFor: valueOptions,
    onEdit: vi.fn(),
    onAdd: vi.fn(),
    onRemove: vi.fn(),
    onClose: vi.fn(),
    ...over,
  };
}

beforeEach(() => {
  vi.clearAllMocks();
});

describe("<ManageTagsModal>", () => {
  // ── Dialog semantics ─────────────────────────────────────────────────────────

  it("renders with role=dialog and aria-modal=true", () => {
    render(<ManageTagsModal {...baseProps()} />);
    const dialog = screen.getByRole("dialog");
    expect(dialog).toBeInTheDocument();
    expect(dialog).toHaveAttribute("aria-modal", "true");
  });

  it("is labelled by the title (sample name) via aria-labelledby", () => {
    render(<ManageTagsModal {...baseProps()} />);
    const dialog = screen.getByRole("dialog");
    const labelledById = dialog.getAttribute("aria-labelledby");
    expect(labelledById).toBeTruthy();
    const title = document.getElementById(labelledById!);
    expect(title).toHaveTextContent("JC042 — LL37");
  });

  it("renders nothing when open=false", () => {
    render(<ManageTagsModal {...baseProps({ open: false })} />);
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  // ── Head + Footer ─────────────────────────────────────────────────────────────

  it("shows the kicker 'Sample tags' and the sample name as h2", () => {
    render(<ManageTagsModal {...baseProps()} />);
    expect(screen.getByText("Sample tags")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "JC042 — LL37" })).toBeInTheDocument();
  });

  it("shows the footer note and a Done button", () => {
    render(<ManageTagsModal {...baseProps()} />);
    expect(screen.getByText(/Tags apply to the whole sample/)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Done" })).toBeInTheDocument();
  });

  it("Done button calls onClose", () => {
    const onClose = vi.fn();
    render(<ManageTagsModal {...baseProps({ onClose })} />);
    fireEvent.click(screen.getByRole("button", { name: "Done" }));
    expect(onClose).toHaveBeenCalledOnce();
  });

  // ── Tag rows ──────────────────────────────────────────────────────────────────

  it("renders one row per tag", () => {
    render(<ManageTagsModal {...baseProps()} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    expect(rows).toHaveLength(TAGS.length);
  });

  it("each row has a key TagSuggest, value TagSuggest, and remove button", () => {
    render(<ManageTagsModal {...baseProps()} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    const firstRow = rows[0]!;
    const suggests = within(firstRow).getAllByTestId("tag-suggest");
    expect(suggests).toHaveLength(2);
    expect(suggests[0]).toHaveAttribute("data-mode", "key");
    expect(suggests[1]).toHaveAttribute("data-mode", "value");
    // Remove button with composed aria-label
    expect(within(firstRow).getByRole("button", { name: /remove lipid LL37/i })).toBeInTheDocument();
  });

  // ── Edit commits by id ────────────────────────────────────────────────────────

  it("editing a value fires onEdit(tagId, key, newValue)", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    const firstRow = rows[0]!;
    const suggests = within(firstRow).getAllByTestId("tag-suggest");
    // The value combobox (mode=value)
    const valueCombobox = within(suggests[1]!).getByRole("combobox");
    // Type into it and press Enter to commit
    fireEvent.change(valueCombobox, { target: { value: "DMPC" } });
    fireEvent.keyDown(valueCombobox, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith(1, "lipid", "DMPC");
  });

  it("editing a key fires onEdit(tagId, newKey, value)", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    const secondRow = rows[1]!;
    const suggests = within(secondRow).getAllByTestId("tag-suggest");
    const keyCombobox = within(suggests[0]!).getByRole("combobox");
    fireEvent.change(keyCombobox, { target: { value: "dose" } });
    fireEvent.keyDown(keyCombobox, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith(2, "dose", "37C");
  });

  // ── Delete by id ──────────────────────────────────────────────────────────────

  it("remove button calls onRemove with the correct tag id", () => {
    const onRemove = vi.fn();
    render(<ManageTagsModal {...baseProps({ onRemove })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    fireEvent.click(within(rows[0]!).getByRole("button", { name: /remove lipid LL37/i }));
    expect(onRemove).toHaveBeenCalledWith(1);
    expect(onRemove).not.toHaveBeenCalledWith(2);
  });

  it("clicking the second remove button calls onRemove with id 2", () => {
    const onRemove = vi.fn();
    render(<ManageTagsModal {...baseProps({ onRemove })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    fireEvent.click(within(rows[1]!).getByRole("button", { name: /remove temperature 37C/i }));
    expect(onRemove).toHaveBeenCalledWith(2);
  });

  // ── Add row: existing key rejection ──────────────────────────────────────────

  it("adding a key that already exists shows a role=alert, does not call onAdd", () => {
    const onAdd = vi.fn();
    render(<ManageTagsModal {...baseProps({ onAdd })} />);
    const addRow = screen.getByTestId("manage-tag-add-row");
    const addKeyCombobox = within(addRow).getAllByRole("combobox")[0]!;
    // Type "lipid" which already exists in TAGS
    fireEvent.change(addKeyCombobox, { target: { value: "lipid" } });
    fireEvent.keyDown(addKeyCombobox, { key: "Enter" });
    // Then click "Add"
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(onAdd).not.toHaveBeenCalled();
    const alert = screen.getByRole("alert");
    expect(alert).toHaveTextContent(/already has a tag with that key/);
  });

  it("the duplicate-key alert is cleared when the key is changed", () => {
    const onAdd = vi.fn();
    render(<ManageTagsModal {...baseProps({ onAdd })} />);
    const addRow = screen.getByTestId("manage-tag-add-row");
    const addKeyCombobox = within(addRow).getAllByRole("combobox")[0]!;
    // Trigger dup error
    fireEvent.change(addKeyCombobox, { target: { value: "lipid" } });
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(screen.getByRole("alert")).toBeInTheDocument();
    // Fix it by typing a new key
    fireEvent.change(addKeyCombobox, { target: { value: "newkey" } });
    expect(screen.queryByRole("alert")).toBeNull();
  });

  it("adds a new tag when key is unique", () => {
    const onAdd = vi.fn();
    render(<ManageTagsModal {...baseProps({ onAdd })} />);
    const addRow = screen.getByTestId("manage-tag-add-row");
    const [addKeyCombo, addValCombo] = within(addRow).getAllByRole("combobox");
    fireEvent.change(addKeyCombo!, { target: { value: "dose" } });
    fireEvent.change(addValCombo!, { target: { value: "10" } });
    fireEvent.click(screen.getByRole("button", { name: "Add" }));
    expect(onAdd).toHaveBeenCalledWith("dose", "10");
  });

  // ── Focus restore on close ────────────────────────────────────────────────────

  it("restores focus to triggerRef.current when closed", () => {
    const triggerEl = document.createElement("button");
    document.body.appendChild(triggerEl);
    const triggerRef = createRef<HTMLButtonElement>();
    (triggerRef as React.MutableRefObject<HTMLButtonElement>).current = triggerEl;
    triggerEl.focus = vi.fn();

    const { rerender } = render(
      <ManageTagsModal {...baseProps({ triggerRef, open: true })} />,
    );

    act(() => {
      rerender(<ManageTagsModal {...baseProps({ triggerRef, open: false })} />);
    });

    expect(triggerEl.focus).toHaveBeenCalled();
    triggerEl.remove();
  });

  it("does NOT steal focus to the trigger on initial mount with open=false", () => {
    // The modal mounts on every loupe load with open=false; the focus-restore
    // effect must not fire on that first render (no true→false transition).
    const triggerEl = document.createElement("button");
    document.body.appendChild(triggerEl);
    const triggerRef = createRef<HTMLButtonElement>();
    (triggerRef as React.MutableRefObject<HTMLButtonElement>).current = triggerEl;
    triggerEl.focus = vi.fn();

    render(<ManageTagsModal {...baseProps({ triggerRef, open: false })} />);

    expect(triggerEl.focus).not.toHaveBeenCalled();
    triggerEl.remove();
  });

  // ── Edit collision: single-valued-key rule on the edit path ──────────────────

  it("editing a key to collide with ANOTHER row's key is rejected inline (role=alert), no onEdit", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    // TAGS = [{id:1, lipid/LL37}, {id:2, temperature/37C}]. Rename row 2's key
    // to "lipid" — collides with row 1.
    const rows = screen.getAllByTestId("manage-tag-row");
    const secondRow = rows[1]!;
    const keyCombobox = within(within(secondRow).getAllByTestId("tag-suggest")[0]!).getByRole("combobox");
    fireEvent.change(keyCombobox, { target: { value: "lipid" } });
    fireEvent.keyDown(keyCombobox, { key: "Enter" });
    expect(onEdit).not.toHaveBeenCalled();
    const alert = screen.getByTestId("manage-tag-edit-dup-error");
    expect(alert).toHaveAttribute("role", "alert");
    expect(alert).toHaveTextContent(/already has a tag with that key/);
    // The colliding key field is marked invalid + tied to the error.
    expect(keyCombobox).toHaveAttribute("aria-invalid", "true");
    expect(keyCombobox).toHaveAttribute("aria-describedby", alert.getAttribute("id"));
  });

  it("the edit-collision alert clears when the key is changed again", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    const keyCombobox = within(within(rows[1]!).getAllByTestId("tag-suggest")[0]!).getByRole("combobox");
    fireEvent.change(keyCombobox, { target: { value: "lipid" } });
    fireEvent.keyDown(keyCombobox, { key: "Enter" });
    expect(screen.getByTestId("manage-tag-edit-dup-error")).toBeInTheDocument();
    // Typing a different key clears the per-row error.
    fireEvent.change(keyCombobox, { target: { value: "dose" } });
    expect(screen.queryByTestId("manage-tag-edit-dup-error")).toBeNull();
  });

  it("editing a key to a NON-colliding value still commits via onEdit", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    const keyCombobox = within(within(rows[1]!).getAllByTestId("tag-suggest")[0]!).getByRole("combobox");
    fireEvent.change(keyCombobox, { target: { value: "dose" } });
    fireEvent.keyDown(keyCombobox, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith(2, "dose", "37C");
    expect(screen.queryByTestId("manage-tag-edit-dup-error")).toBeNull();
  });

  it("re-committing a row's OWN unchanged key is allowed (self-collision excluded)", () => {
    const onEdit = vi.fn();
    render(<ManageTagsModal {...baseProps({ onEdit })} />);
    const rows = screen.getAllByTestId("manage-tag-row");
    // Row 1's key is "lipid"; committing it unchanged must not flag a collision.
    const keyCombobox = within(within(rows[0]!).getAllByTestId("tag-suggest")[0]!).getByRole("combobox");
    fireEvent.keyDown(keyCombobox, { key: "Enter" });
    expect(screen.queryByTestId("manage-tag-edit-dup-error")).toBeNull();
    expect(onEdit).toHaveBeenCalledWith(1, "lipid", "LL37");
  });
});
