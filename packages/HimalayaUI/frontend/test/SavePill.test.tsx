import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SavePill } from "../src/components/SavePill";

describe("SavePill — Compare UX C-7", () => {
  it("hides when not dirty", () => {
    const { queryByTestId } = render(
      <SavePill dirty={false} mode={{ kind: "viewing" }} onSave={() => {}} isSaving={false}/>);
    expect(queryByTestId("save-pill")).toBeNull();
  });

  it("shows 'Save changes' for editing-mine", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Save changes");
  });

  it("shows 'Save as fork…' for editing-as-fork-of", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-as-fork-of", parentId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent(/Save as fork/);
  });

  it("shows 'Save' for creating-blank", () => {
    render(<SavePill dirty={true} mode={{ kind: "creating-blank" }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Save");
  });

  it("shows 'Save fork' for creating-from-fork (post-morph)", () => {
    render(<SavePill dirty={true} mode={{ kind: "creating-from-fork", parentId: 1 }} onSave={() => {}} isSaving={false}/>);
    expect(screen.getByTestId("save-pill")).toHaveTextContent(/Save fork/);
  });

  it("calls onSave when clicked", () => {
    const onSave = vi.fn();
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={onSave} isSaving={false}/>);
    fireEvent.click(screen.getByTestId("save-pill"));
    expect(onSave).toHaveBeenCalled();
  });

  it("disables and shows 'Saving…' while isSaving", () => {
    render(<SavePill dirty={true} mode={{ kind: "editing-mine", draftId: 1 }} onSave={() => {}} isSaving={true}/>);
    expect(screen.getByTestId("save-pill")).toBeDisabled();
    expect(screen.getByTestId("save-pill")).toHaveTextContent("Saving…");
  });
});
