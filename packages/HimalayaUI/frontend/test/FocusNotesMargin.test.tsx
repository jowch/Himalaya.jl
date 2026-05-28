import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { Sample } from "../src/api";
import { FocusNotesMargin } from "../src/components/FocusNotesMargin";

const SAMPLE: Sample = {
  id: 1, experiment_id: 1, name: "JC001",
  display_name: "Sample 1", notes: "weak shoulder at q≈0.064", tags: [],
};

beforeEach(() => {
  localStorage.clear();
});

describe("FocusNotesMargin", () => {
  it("renders the sample's notes", () => {
    render(<FocusNotesMargin sample={SAMPLE} onSaveNotes={() => {}} />);
    expect(screen.getByTestId("focus-notes-margin")).toBeInTheDocument();
    expect(screen.getByDisplayValue("weak shoulder at q≈0.064")).toBeInTheDocument();
  });

  it("does not clobber a mid-edit draft when the sample prop updates", async () => {
    const user = userEvent.setup();
    const { rerender } = render(
      <FocusNotesMargin sample={SAMPLE} onSaveNotes={() => {}} />,
    );
    const ta = screen.getByTestId("focus-notes-input") as HTMLTextAreaElement;
    await user.click(ta);
    await user.clear(ta);
    await user.type(ta, "my in-progress edit");
    // a background refetch delivers a different notes value while focused:
    rerender(
      <FocusNotesMargin
        sample={{ ...SAMPLE, notes: "someone else's save" }}
        onSaveNotes={() => {}}
      />,
    );
    expect(ta.value).toBe("my in-progress edit"); // focus-gate held
  });

  it("R3-F04: shows a count badge when notes are present", () => {
    render(<FocusNotesMargin sample={SAMPLE} onSaveNotes={() => {}} />);
    expect(screen.getByTestId("focus-notes-count")).toHaveTextContent("1");
  });

  it("R3-F04: hides the count badge when notes are empty", () => {
    render(<FocusNotesMargin sample={{ ...SAMPLE, notes: "" }} onSaveNotes={() => {}} />);
    expect(screen.queryByTestId("focus-notes-count")).toBeNull();
  });

  it("R3-F04: mono-formats q ≈ N.NNN references in the note body", () => {
    render(
      <FocusNotesMargin
        sample={{ ...SAMPLE, notes: "weak shoulder at q ≈ 0.064 here" }}
        onSaveNotes={() => {}}
      />,
    );
    const refs = screen.getAllByTestId("focus-notes-qref");
    expect(refs.length).toBeGreaterThanOrEqual(1);
    expect(refs[0]).toHaveTextContent("q ≈ 0.064");
  });

  it("R3-F04: renders the dashed add-a-note affordance", () => {
    render(<FocusNotesMargin sample={{ ...SAMPLE, notes: "" }} onSaveNotes={() => {}} />);
    expect(screen.getByTestId("focus-notes-add")).toBeInTheDocument();
  });

  it("calls onSaveNotes with the edited value on blur", async () => {
    const user = userEvent.setup();
    const onSaveNotes = vi.fn();
    render(<FocusNotesMargin sample={SAMPLE} onSaveNotes={onSaveNotes} />);
    const ta = screen.getByTestId("focus-notes-input");
    await user.click(ta);
    await user.clear(ta);
    await user.type(ta, "updated note");
    await user.tab(); // blur
    expect(onSaveNotes).toHaveBeenCalledWith("updated note");
  });
});
