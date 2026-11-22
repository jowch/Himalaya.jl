// test/DirectoryPickerField.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { DirectoryPickerField } from "../src/components/DirectoryPickerField";

describe("DirectoryPickerField (Phase E1)", () => {
  it("shows suggestions and Tab completes the top one", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/Volumes/data/ssrl/2026_04/1p7m", "/Volumes/data/ssrl/2026_04/10x"]}
      />,
    );
    // Input puts data-testid on its WRAPPER; onKeyDown (via ...rest) lands on
    // the inner <input>. Fire on the inner input so the handler runs.
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Tab" });
    expect(onChange).toHaveBeenCalledWith("/Volumes/data/ssrl/2026_04/1p7m/");
  });

  it("arrow-down moves the active suggestion, Enter completes it", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/a/one", "/a/two"]}
      />,
    );
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "ArrowDown" });
    fireEvent.keyDown(input, { key: "ArrowDown" });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith("/a/two/");
  });

  it("highlights no suggestion until the user moves into the list", () => {
    render(<DirectoryPickerField value="/a/o" onChange={() => {}} suggestions={["/a/one", "/a/two"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input);
    expect(screen.getAllByRole("option").filter((o) => o.getAttribute("aria-selected") === "true")).toHaveLength(0);
  });

  it("Enter with no suggestion highlighted submits the typed folder instead of completing a subfolder", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/beamtime/setup1" onChange={onChange} onSubmit={onSubmit} suggestions={["/beamtime/setup1/analysis", "/beamtime/setup1/data"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onSubmit).toHaveBeenCalledTimes(1);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("Enter submits when there are no suggestions at all", () => {
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/beamtime/setup1/data" onChange={() => {}} onSubmit={onSubmit} suggestions={[]} />);
    fireEvent.keyDown(screen.getByTestId("dirpicker-input").querySelector("input")!, { key: "Enter" });
    expect(onSubmit).toHaveBeenCalledTimes(1);
  });

  it("opening the suggestion popover keeps focus in the input, not the dialog panel", () => {
    // Regression: Popover's open-effect focuses its own panel by default. Without
    // initialFocusRef → the input, clicking the field stole focus onto the
    // role=dialog panel, so typing and Tab-completion never reached the input's
    // onKeyDown. (The other tests fire keyDown ON the input directly, so they
    // pass even when focus is actually on the panel — they cannot catch this.)
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={() => {}}
        suggestions={["/Volumes/data/ssrl/2026_04/1p7m"]}
      />,
    );
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    // Click the field — this is the trigger that opens the Popover.
    fireEvent.click(input);
    // The suggestion dialog is open …
    expect(screen.getByRole("dialog")).toBeInTheDocument();
    // … and focus is on the INPUT (so the caret stays, typing works, and the
    // onKeyDown Tab/↑↓ handlers fire), NOT on the panel.
    expect(document.activeElement).toBe(input);
  });

  it("shows each suggestion's basename but completes the full path", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Users/me/projects/himalaya-devdata/mi"
        onChange={onChange}
        suggestions={["/Users/me/projects/himalaya-devdata/mini-1p7m"]}
      />,
    );
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input); // open the popover so the list renders
    const opt = screen.getByRole("option");
    // Label is the basename, not the full path …
    expect(opt).toHaveTextContent("mini-1p7m");
    expect(opt.textContent).not.toContain("/");
    // … the full path stays available as a title for disambiguation …
    expect(opt).toHaveAttribute("title", "/Users/me/projects/himalaya-devdata/mini-1p7m");
    // … and selecting it completes the FULL path into the field.
    fireEvent.click(opt);
    expect(onChange).toHaveBeenCalledWith("/Users/me/projects/himalaya-devdata/mini-1p7m/");
  });

  it("Shift+Enter with nothing highlighted neither completes nor submits", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/a/o" onChange={onChange} onSubmit={onSubmit} suggestions={["/a/one", "/a/two"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    const errors: unknown[] = [];
    const onError = (e: ErrorEvent): void => {
      errors.push(e.error);
      e.preventDefault();
    };
    window.addEventListener("error", onError);
    try {
      fireEvent.keyDown(input, { key: "Enter", shiftKey: true });
      fireEvent.keyDown(input, { key: "Enter", altKey: true });
    } catch (e) {
      errors.push(e);
    } finally {
      window.removeEventListener("error", onError);
    }
    expect(errors).toEqual([]);
    expect(onChange).not.toHaveBeenCalled();
    expect(onSubmit).not.toHaveBeenCalled();
  });

  it("Tab leaves the field when the typed folder is complete and the list only holds its subfolders", () => {
    const onChange = vi.fn();
    render(<DirectoryPickerField value="/b/day1/setup2" onChange={onChange} suggestions={["/b/day1/setup2/analysis", "/b/day1/setup2/data"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    expect(fireEvent.keyDown(input, { key: "Tab" })).toBe(true);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("Enter that confirms an IME composition does not submit", () => {
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/a/o" onChange={() => {}} onSubmit={onSubmit} suggestions={[]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Enter", isComposing: true });
    expect(onSubmit).not.toHaveBeenCalled();
  });
  it("Tab does not complete a suggestion left over from an earlier prefix", () => {
    const onChange = vi.fn();
    render(<DirectoryPickerField value="/b/day2" onChange={onChange} suggestions={["/b/day1", "/b/day1x"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    expect(fireEvent.keyDown(input, { key: "Tab" })).toBe(true);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("Tab skips leftover suggestions and completes the first that matches what is typed", () => {
    const onChange = vi.fn();
    render(<DirectoryPickerField value="/b/se" onChange={onChange} suggestions={["/b/notes", "/b/setup1"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Tab" });
    expect(onChange).toHaveBeenCalledWith("/b/setup1/");
  });

  it("keeps the suggestions out of the Tab order", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} suggestions={["/b/setup1", "/b/setup2"]} />);
    fireEvent.click(screen.getByTestId("dirpicker-input").querySelector("input")!);
    const options = screen.getAllByRole("option");
    expect(options).toHaveLength(2);
    for (const o of options) expect(o).toHaveAttribute("tabindex", "-1");
  });
  it("points the field at the highlighted suggestion for screen readers", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} suggestions={["/b/setup1", "/b/setup2"]} />);
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.click(input);
    expect(input).not.toHaveAttribute("aria-activedescendant");
    fireEvent.keyDown(input, { key: "ArrowDown" });
    fireEvent.keyDown(input, { key: "ArrowDown" });
    const id = input.getAttribute("aria-activedescendant")!;
    expect(document.getElementById(id)).toBe(screen.getAllByRole("option")[1]);
  });
});
