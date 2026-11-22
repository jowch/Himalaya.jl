import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { DirectoryPickerField, folderHint } from "../src/components/DirectoryPickerField";
import type { FolderSuggestion } from "../src/api";

const f = (path: string, extra: Partial<FolderSuggestion> = {}): FolderSuggestion =>
  ({ path, kind: null, setups: null, added: null, described: true, ...extra });

const input = (): HTMLInputElement => screen.getByRole("combobox", { name: "Data directory" });
const selected = (): string[] =>
  screen.queryAllByRole("option").filter((o) => o.getAttribute("aria-selected") === "true").map((o) => o.textContent ?? "");

describe("DirectoryPickerField", () => {
  it("shows the folders as soon as the field has focus, without a click", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1"), f("/b/setup2")]} autoFocus />);
    expect(document.activeElement).toBe(input());
    expect(screen.getAllByRole("option").map((o) => o.textContent)).toEqual(["setup1", "setup2"]);
    expect(input()).toHaveAttribute("aria-expanded", "true");
  });

  it("hides the list when focus leaves the field and shows it again on return", () => {
    render(<><DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus /><button>next</button></>);
    fireEvent.blur(input(), { relatedTarget: screen.getByRole("button", { name: "next" }) });
    expect(screen.queryByRole("listbox")).toBeNull();
    fireEvent.focus(input());
    expect(screen.getByRole("listbox")).toBeInTheDocument();
  });

  it("highlights nothing until Tab or ↓ moves into the list", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1"), f("/b/setup2")]} autoFocus />);
    expect(selected()).toEqual([]);
    expect(input()).not.toHaveAttribute("aria-activedescendant");
    fireEvent.keyDown(input(), { key: "Tab" });
    expect(selected()).toEqual(["setup1"]);
    fireEvent.keyDown(input(), { key: "ArrowDown" });
    expect(selected()).toEqual(["setup2"]);
    expect(document.getElementById(input().getAttribute("aria-activedescendant")!)).toBe(screen.getAllByRole("option")[1]);
  });

  it("Tab and Shift+Tab move through the list and wrap", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1"), f("/b/setup2")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Tab" });
    expect(selected()).toEqual(["setup1"]);
    fireEvent.keyDown(input(), { key: "Tab", shiftKey: true });
    expect(selected()).toEqual(["setup2"]);
  });

  it("Enter on a highlighted folder fills it in with a slash, to go into it", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/b/s" onChange={onChange} onSubmit={onSubmit} folders={[f("/b/setup1"), f("/b/setup2")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith("/b/setup2/");
    expect(onSubmit).not.toHaveBeenCalled();
  });

  it("Enter with nothing highlighted submits the typed folder", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/b/day1/" onChange={onChange} onSubmit={onSubmit} folders={[f("/b/day1/setup1")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onSubmit).toHaveBeenCalledTimes(1);
    expect(onChange).not.toHaveBeenCalled();
  });

  it("Enter submits when there are no folders at all", () => {
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/b/day1" onChange={() => {}} onSubmit={onSubmit} folders={[]} />);
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onSubmit).toHaveBeenCalledTimes(1);
  });

  it("Shift+Enter, Alt+Enter and an IME Enter do nothing", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/a/o" onChange={onChange} onSubmit={onSubmit} folders={[f("/a/one")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Enter", shiftKey: true });
    fireEvent.keyDown(input(), { key: "Enter", altKey: true });
    fireEvent.keyDown(input(), { key: "Enter", isComposing: true });
    expect(onChange).not.toHaveBeenCalled();
    expect(onSubmit).not.toHaveBeenCalled();
  });

  it("⌘/Ctrl+Enter is left to the page, even with a folder highlighted", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/a/o" onChange={onChange} onSubmit={onSubmit} folders={[f("/a/one")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Tab" });
    expect(fireEvent.keyDown(input(), { key: "Enter", metaKey: true })).toBe(true);
    expect(fireEvent.keyDown(input(), { key: "Enter", ctrlKey: true })).toBe(true);
    expect(onChange).not.toHaveBeenCalled();
    expect(onSubmit).not.toHaveBeenCalled();
  });

  it("Esc closes the list in one press, and the next Esc reaches the page", () => {
    const onPageKey = vi.fn();
    render(<div onKeyDown={(e) => { if (e.key === "Escape") onPageKey(); }}><DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus /></div>);
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Escape" });
    expect(screen.queryByRole("listbox")).toBeNull();
    expect(onPageKey).not.toHaveBeenCalled();
    fireEvent.keyDown(input(), { key: "Escape" });
    expect(onPageKey).toHaveBeenCalledTimes(1);
  });

  it("typing clears the highlight, so Enter adds the typed folder rather than a stale row", () => {
    const onChange = vi.fn();
    const onSubmit = vi.fn();
    render(<DirectoryPickerField value="/b/be" onChange={onChange} onSubmit={onSubmit} folders={[f("/b/beamtime-10"), f("/b/beamtime-11")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.change(input(), { target: { value: "/b/bea" } });
    expect(selected()).toEqual([]);
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onSubmit).toHaveBeenCalledTimes(1);
    expect(onChange).toHaveBeenCalledTimes(1);
  });

  it("after a folder is filled in, its subfolders stay hidden until asked for, so Tab leaves the field", () => {
    const onChange = vi.fn();
    const { rerender } = render(<DirectoryPickerField value="/b/be" onChange={onChange} folders={[f("/b/beamtime")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Tab" });
    fireEvent.keyDown(input(), { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith("/b/beamtime/");
    rerender(<DirectoryPickerField value="/b/beamtime/" onChange={onChange} folders={[f("/b/beamtime/setup1"), f("/b/beamtime/setup2")]} autoFocus />);
    expect(screen.queryByRole("listbox")).toBeNull();
    expect(fireEvent.keyDown(input(), { key: "Tab" })).toBe(true);
    fireEvent.keyDown(input(), { key: "ArrowDown" });
    expect(selected()).toEqual(["setup1"]);
  });

  it("a click on a folder fills it in and closes the list", () => {
    const onChange = vi.fn();
    render(<DirectoryPickerField value="/b/" onChange={onChange} folders={[f("/b/one")]} autoFocus />);
    fireEvent.click(screen.getByRole("option"));
    expect(onChange).toHaveBeenCalledWith("/b/one/");
    expect(screen.queryByRole("listbox")).toBeNull();
  });

  it("marks a folder Himalaya ran out of time to check", () => {
    render(<DirectoryPickerField value="/b/" onChange={() => {}} folders={[f("/b/slow", { described: false })]} autoFocus />);
    expect(screen.getByRole("option")).toHaveAttribute("aria-label", "slow, not checked");
  });

  it("after Esc, Tab leaves the field and ↓ opens the list again", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Escape" });
    expect(fireEvent.keyDown(input(), { key: "Tab" })).toBe(true);
    fireEvent.keyDown(input(), { key: "ArrowDown" });
    expect(selected()).toEqual(["setup1"]);
  });

  it("typing again after Esc shows the list again", () => {
    const { rerender } = render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus />);
    fireEvent.keyDown(input(), { key: "Escape" });
    fireEvent.change(input(), { target: { value: "/b/se" } });
    rerender(<DirectoryPickerField value="/b/se" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus />);
    expect(screen.getByRole("listbox")).toBeInTheDocument();
  });

  it("Shift+Tab with nothing highlighted leaves the field", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1")]} autoFocus />);
    expect(fireEvent.keyDown(input(), { key: "Tab", shiftKey: true })).toBe(true);
  });

  it("hides a list that only repeats the folder already typed, so Tab leaves the field", () => {
    render(<DirectoryPickerField value="/b/day1" onChange={() => {}} folders={[f("/b/day1", { kind: "beamtime", setups: 2 })]} autoFocus />);
    expect(screen.queryByRole("listbox")).toBeNull();
    expect(fireEvent.keyDown(input(), { key: "Tab" })).toBe(true);
  });

  it("shows the folder name and what it holds, and a click fills in the full path", async () => {
    const onChange = vi.fn();
    render(<DirectoryPickerField value="/b/" onChange={onChange} folders={[
      f("/b/2026-11", { kind: "beamtime", setups: 2 }),
      f("/b/setup1", { kind: "setup", added: { id: 3, name: "Oct" } }),
      f("/b/notes"),
    ]} autoFocus />);
    const options = screen.getAllByRole("option");
    expect(options.map((o) => o.getAttribute("aria-label"))).toEqual([
      "2026-11, beamtime, 2 setups",
      "setup1, setup · added as Oct",
      "notes",
    ]);
    fireEvent.mouseEnter(options[0]!);
    expect(await screen.findByRole("tooltip")).toHaveTextContent("/b/2026-11");
    fireEvent.click(options[0]!);
    expect(onChange).toHaveBeenCalledWith("/b/2026-11/");
  });

  it("says when more folders match than are listed", () => {
    render(<DirectoryPickerField value="/b/" onChange={() => {}} folders={[f("/b/a")]} more autoFocus />);
    expect(screen.getByTestId("dirpicker-suggestions-note")).toHaveTextContent("Keep typing");
  });

  it("keeps the options out of the Tab order and focus in the field when one is pressed", () => {
    render(<DirectoryPickerField value="/b/s" onChange={() => {}} folders={[f("/b/setup1"), f("/b/setup2")]} autoFocus />);
    for (const o of screen.getAllByRole("option")) expect(o).toHaveAttribute("tabindex", "-1");
    expect(fireEvent.mouseDown(screen.getAllByRole("option")[0]!)).toBe(false);
  });

  it("folderHint words a beamtime with one setup in the singular", () => {
    expect(folderHint(f("/b/x", { kind: "beamtime", setups: 1 }))).toBe("beamtime, 1 setup");
    expect(folderHint(f("/b/x"))).toBeUndefined();
  });
});
