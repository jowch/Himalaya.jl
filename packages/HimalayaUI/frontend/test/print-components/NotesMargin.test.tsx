import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { NotesMargin } from "../../src/print/components/NotesMargin";

describe("NotesMargin", () => {
  // 1. notes={null} → no body, no count badge; textarea + "Notes" head present
  it("renders head and textarea but no body or badge when notes is null", () => {
    render(<NotesMargin notes={null} onSaveNotes={() => {}} />);
    expect(screen.getByTestId("focus-notes-margin")).toBeInTheDocument();
    // "Notes" kicker text present
    expect(screen.getByText("Notes")).toBeInTheDocument();
    // textarea present
    expect(screen.getByTestId("focus-notes-input")).toBeInTheDocument();
    // no body
    expect(screen.queryByTestId("focus-notes-body")).toBeNull();
    // no count badge
    expect(screen.queryByTestId("focus-notes-count")).toBeNull();
  });

  // 2. notes with two q-refs → body present, count badge "1", two qref spans
  it("renders body with two qref spans and count badge when notes has q-refs", () => {
    render(
      <NotesMargin
        notes="see q ≈ 0.064 and q≈0.128"
        onSaveNotes={() => {}}
      />,
    );
    expect(screen.getByTestId("focus-notes-body")).toBeInTheDocument();
    const badge = screen.getByTestId("focus-notes-count");
    expect(badge.textContent).toBe("1");
    const qrefs = screen.getAllByTestId("focus-notes-qref");
    expect(qrefs).toHaveLength(2);
    expect(qrefs[0].textContent).toContain("q");
    expect(qrefs[1].textContent).toContain("q");
  });

  // 3. typing then blurring fires onSaveNotes with the typed draft
  it("fires onSaveNotes with the typed draft on blur", () => {
    const onSave = vi.fn();
    render(<NotesMargin notes={null} onSaveNotes={onSave} />);
    const ta = screen.getByTestId("focus-notes-input");
    fireEvent.change(ta, { target: { value: "new note text" } });
    fireEvent.blur(ta);
    expect(onSave).toHaveBeenCalledOnce();
    expect(onSave).toHaveBeenCalledWith("new note text");
  });

  // 4. focus-gate: draft is NOT clobbered while focused; re-syncs after blur
  it("does not clobber draft while textarea is focused, re-syncs after blur", () => {
    const onSave = vi.fn();
    const { rerender } = render(
      <NotesMargin notes="original" onSaveNotes={onSave} />,
    );
    const ta = screen.getByTestId("focus-notes-input") as HTMLTextAreaElement;

    // Type something and focus
    fireEvent.focus(ta);
    fireEvent.change(ta, { target: { value: "mid-edit draft" } });

    // Rerender with new notes while still focused — draft must NOT be clobbered
    rerender(<NotesMargin notes="server update" onSaveNotes={onSave} />);
    expect(ta.value).toBe("mid-edit draft");

    // After blur, rerender with new notes should re-sync
    fireEvent.blur(ta);
    rerender(<NotesMargin notes="server update" onSaveNotes={onSave} />);
    expect(ta.value).toBe("server update");
  });

  // 5. onHoverQ: mouseEnter fires with parsed q, mouseLeave fires with undefined
  it("fires onHoverQ(q) on mouseEnter and onHoverQ(undefined) on mouseLeave", () => {
    const onHoverQ = vi.fn();
    render(
      <NotesMargin
        notes="see q ≈ 0.064 for details"
        onSaveNotes={() => {}}
        onHoverQ={onHoverQ}
      />,
    );
    const qref = screen.getByTestId("focus-notes-qref");
    fireEvent.mouseEnter(qref);
    expect(onHoverQ).toHaveBeenCalledWith(0.064);
    fireEvent.mouseLeave(qref);
    expect(onHoverQ).toHaveBeenCalledWith(undefined);
  });
});
