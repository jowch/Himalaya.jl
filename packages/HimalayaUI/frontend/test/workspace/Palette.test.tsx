import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { Palette } from "../../src/workspace/ui/palette/Palette";
import { stepHighlight } from "../../src/workspace/ui/useListHighlight";
import { DEEP_TITLE, cycleOption, rootRows, type PaletteForm, type PaletteGroup, type PalettePage } from "../../src/workspace/ui/palette/model";

const ran: string[] = [];
const run = (id: string) => () => ran.push(id);

function setup(): { groups: PaletteGroup[]; pages: Record<string, PalettePage> } {
  const pages: Record<string, PalettePage> = {
    compare: {
      title: "Compare with",
      arg: "Sample",
      deep: true,
      items: ["B1", "B2", "B3"].map((n) => ({ id: `compare-${n}`, label: n, run: run(`compare ${n}`) })),
    },
    theme: { title: "Theme", arg: "Theme", items: [{ id: "dark", label: "Dark", run: run("dark") }] },
    label: {
      title: "Label",
      arg: "Name",
      items: [{ id: "retake", label: "retake", run: run("retake") }],
      extra: (q) => (q ? [{ id: "create", label: `Create “${q}”`, form: () => newLabel(q) }] : []),
    },
  };
  const groups: PaletteGroup[] = [
    {
      title: "This sample",
      items: [
        { id: "next", label: "Next unreviewed sample", shortcut: "Z", run: run("next") },
        { id: "compare", label: "Compare with a sample", page: "compare" },
      ],
    },
    {
      title: "Experiment",
      items: [
        { id: "home", label: "Back to all experiments", run: run("home") },
        { id: "remove", label: "Remove this experiment", danger: true, typeOnly: true, run: run("remove") },
      ],
    },
    { title: "App", items: [{ id: "theme", label: "Theme", page: "theme" }, { id: "label", label: "Label", page: "label" }] },
    { title: "Samples", searchOnly: true, items: [{ id: "go-b2", label: "Go to B2", run: run("go B2") }] },
  ];
  return { groups, pages };
}

const submitted: Record<string, string>[] = [];

function newLabel(name: string): PaletteForm {
  return {
    title: "New label",
    verb: "create",
    fields: [
      { id: "name", label: "Name", kind: "text", value: name, check: (v) => (v.trim() ? undefined : "Type a name.") },
      { id: "color", label: "Colour", kind: "choice", value: "teal", options: ["rose", "teal", "stone"].map((c) => ({ value: c, label: c[0]!.toUpperCase() + c.slice(1) })) },
    ],
    submit: (v) => submitted.push({ ...v }),
  };
}

function heads(rows: ReturnType<typeof rootRows>): string[] {
  return rows.flatMap((r) => (r.kind === "head" ? [r.title] : []));
}

describe("rootRows", () => {
  it("lists every group but search-only ones, and hides type-only items until typed", () => {
    const { groups, pages } = setup();
    const rows = rootRows(groups, pages, "");
    expect(heads(rows)).toEqual(["This sample", "Experiment", "App"]);
    expect(rows.some((r) => r.kind === "item" && r.item.id === "remove")).toBe(false);
    expect(rootRows(groups, pages, "remove").flatMap((r) => (r.kind === "item" ? [r.item.id] : []))).toEqual(["remove"]);
  });

  it("searches inside deep pages and samples, but not inside other pages", () => {
    const { groups, pages } = setup();
    const rows = rootRows(groups, pages, "b2");
    expect(heads(rows)).toEqual([DEEP_TITLE, "Samples"]);
    expect(rows.flatMap((r) => (r.kind === "item" ? [`${r.path ?? ""}|${r.item.label}`] : []))).toEqual(["Compare with|B2", "|Go to B2"]);
    expect(rootRows(groups, pages, "dark")).toEqual([]);
  });

  it("matches a deep item by its page title and its own words together", () => {
    const { groups, pages } = setup();
    const ids = rootRows(groups, pages, "compare b3").flatMap((r) => (r.kind === "item" ? [r.item.id] : []));
    expect(ids).toEqual(["compare-B3"]);
  });
});

function open(): { onClose: ReturnType<typeof vi.fn>; rerender: (pages: Record<string, PalettePage>) => void } {
  ran.length = 0;
  const onClose = vi.fn();
  const { groups, pages } = setup();
  const view = render(<Palette groups={groups} pages={pages} onClose={onClose} />);
  return { onClose, rerender: (next) => view.rerender(<Palette groups={groups} pages={next} onClose={onClose} />) };
}

const field = (): HTMLElement => screen.getByRole("combobox");
const selected = (): string | null => screen.getAllByRole("option").find((o) => o.getAttribute("aria-selected") === "true")?.getAttribute("data-palette-item") ?? null;

describe("Palette", () => {
  it("opens a page in place, with its title and argument in the bar", () => {
    open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    expect(screen.getByRole("dialog", { name: "Compare with" })).toBeInTheDocument();
    expect(screen.getByRole("combobox", { name: "Sample" })).toHaveFocus();
    expect(screen.getByTestId("workspace-palette-bar")).toHaveTextContent("Compare withSample");
    expect(screen.getAllByRole("option").map((o) => o.textContent)).toEqual(["B1", "B2", "B3"]);
  });

  it("Esc goes back one level and keeps what was typed, and closes only from the top", () => {
    const { onClose } = open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.change(field(), { target: { value: "b" } });
    fireEvent.keyDown(field(), { key: "Escape" });
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toHaveValue("compare");
    expect(onClose).not.toHaveBeenCalled();
    fireEvent.keyDown(field(), { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("Backspace in an empty field goes back, and in a filled one just edits", () => {
    const { onClose } = open();
    fireEvent.click(screen.getByRole("option", { name: "Compare with a sample" }));
    fireEvent.change(field(), { target: { value: "b" } });
    fireEvent.keyDown(field(), { key: "Backspace" });
    expect(screen.getByRole("combobox", { name: "Sample" })).toBeInTheDocument();
    fireEvent.change(field(), { target: { value: "" } });
    fireEvent.keyDown(field(), { key: "Backspace" });
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toBeInTheDocument();
    fireEvent.keyDown(field(), { key: "Backspace" });
    expect(onClose).not.toHaveBeenCalled();
  });

  it("Tab, Shift+Tab and the arrows move the highlight and wrap, with focus kept in the field", () => {
    open();
    expect(selected()).toBe("next");
    fireEvent.keyDown(field(), { key: "Tab", shiftKey: true });
    expect(selected()).toBe("label");
    fireEvent.keyDown(field(), { key: "Tab" });
    expect(selected()).toBe("next");
    fireEvent.keyDown(field(), { key: "ArrowUp" });
    expect(selected()).toBe("label");
    fireEvent.keyDown(field(), { key: "ArrowDown" });
    fireEvent.keyDown(field(), { key: "ArrowDown" });
    expect(selected()).toBe("compare");
    expect(field()).toHaveFocus();
  });

  it("Enter on a page's item runs it and closes the palette", () => {
    const { onClose } = open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.change(field(), { target: { value: "b3" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    expect(ran).toEqual(["compare B3"]);
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("runs a deep match straight from the top, naming the page it lives in", () => {
    open();
    fireEvent.change(field(), { target: { value: "b2" } });
    expect(screen.getAllByRole("option")[0]).toHaveAccessibleName("Compare with › B2");
    fireEvent.keyDown(field(), { key: "Enter" });
    expect(ran).toEqual(["compare B2"]);
  });

  it("styles a type-only danger item as danger once typed", () => {
    open();
    fireEvent.change(field(), { target: { value: "remove" } });
    expect(screen.getByRole("option", { name: "Remove this experiment" }).firstElementChild).toHaveClass("text-error");
  });

  it("turns a Create row into an inline form: Tab moves between fields, arrows change a choice, Enter submits", () => {
    const { onClose } = open();
    submitted.length = 0;
    fireEvent.click(screen.getByRole("option", { name: "Label" }));
    fireEvent.change(field(), { target: { value: "dehydrated" } });
    expect(screen.getAllByRole("option").map((o) => o.getAttribute("aria-label"))).toEqual(["Create “dehydrated”"]);
    fireEvent.keyDown(field(), { key: "Enter" });
    expect(screen.getByRole("dialog", { name: "New label" })).toBeInTheDocument();
    const name = screen.getByRole("textbox", { name: "Name" });
    expect(name).toHaveFocus();
    expect(name).toHaveValue("dehydrated");
    fireEvent.keyDown(name, { key: "Tab" });
    const colour = screen.getByRole("combobox", { name: "Colour" });
    expect(colour).toHaveFocus();
    expect(screen.getByRole("listbox", { name: "Colour" })).toBeInTheDocument();
    expect(screen.getByRole("option", { name: "Teal" })).toHaveAttribute("aria-selected", "true");
    fireEvent.keyDown(colour, { key: "ArrowDown" });
    expect(screen.getByRole("option", { name: "Stone" })).toHaveAttribute("aria-selected", "true");
    fireEvent.keyDown(colour, { key: "Tab" });
    expect(name).toHaveFocus();
    fireEvent.keyDown(name, { key: "Enter" });
    expect(submitted).toEqual([{ name: "dehydrated", color: "stone" }]);
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("a held Enter in a form submits once, and a held Esc steps back one level", () => {
    const { onClose } = open();
    submitted.length = 0;
    fireEvent.click(screen.getByRole("option", { name: "Label" }));
    fireEvent.change(field(), { target: { value: "x" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    const name = screen.getByRole("textbox", { name: "Name" });
    fireEvent.keyDown(name, { key: "Enter", repeat: true });
    expect(submitted).toEqual([]);
    fireEvent.keyDown(name, { key: "Escape" });
    fireEvent.keyDown(screen.getByRole("combobox"), { key: "Escape", repeat: true });
    expect(screen.getByRole("combobox", { name: "Name" })).toHaveValue("x");
    expect(onClose).not.toHaveBeenCalled();
  });

  it("refuses to submit an invalid field, outlines it and says why, and Esc returns to the list as it was", () => {
    open();
    submitted.length = 0;
    fireEvent.click(screen.getByRole("option", { name: "Label" }));
    fireEvent.change(field(), { target: { value: "x" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    const name = screen.getByRole("textbox", { name: "Name" });
    fireEvent.change(name, { target: { value: " " } });
    fireEvent.keyDown(name, { key: "Enter" });
    expect(submitted).toEqual([]);
    expect(name).toHaveAttribute("aria-invalid", "true");
    expect(name).toHaveAccessibleDescription("Type a name.");
    fireEvent.keyDown(name, { key: "Escape" });
    expect(screen.getByRole("combobox", { name: "Name" })).toHaveValue("x");
  });

  it("Backspace on a choice or an empty text field steps back a field, then out of the form", () => {
    open();
    fireEvent.click(screen.getByRole("option", { name: "Label" }));
    fireEvent.change(field(), { target: { value: "x" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.keyDown(screen.getByRole("textbox", { name: "Name" }), { key: "Tab" });
    fireEvent.keyDown(screen.getByRole("combobox", { name: "Colour" }), { key: "Backspace" });
    const name = screen.getByRole("textbox", { name: "Name" });
    expect(name).toHaveFocus();
    fireEvent.keyDown(name, { key: "Backspace" });
    expect(screen.getByRole("textbox", { name: "Name" })).toBeInTheDocument();
    fireEvent.change(name, { target: { value: "" } });
    fireEvent.keyDown(name, { key: "Backspace" });
    expect(screen.getByRole("dialog", { name: "Label" })).toBeInTheDocument();
  });

  it("keeps focus in the field when the footer, a group head or the empty space is clicked", () => {
    open();
    const footer = screen.getByText("move");
    expect(fireEvent.mouseDown(footer)).toBe(false);
    expect(fireEvent.mouseDown(screen.getByText("This sample"))).toBe(false);
    expect(fireEvent.mouseDown(field())).toBe(true);
    expect(field()).toHaveFocus();
    fireEvent.keyDown(field(), { key: "Escape" });
  });

  it("ignores a held Enter, so opening a page never also runs its first row", () => {
    const { onClose } = open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.keyDown(field(), { key: "Enter", repeat: true });
    fireEvent.keyDown(field(), { key: "Enter", repeat: true });
    expect(ran).toEqual([]);
    expect(onClose).not.toHaveBeenCalled();
    expect(screen.getByRole("dialog", { name: "Compare with" })).toBeInTheDocument();
  });

  it("ignores a held Esc or Backspace, so one press goes back exactly one level", () => {
    const { onClose } = open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.keyDown(field(), { key: "Escape" });
    fireEvent.keyDown(field(), { key: "Escape", repeat: true });
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toHaveValue("compare");
    fireEvent.keyDown(field(), { key: "Enter" });
    fireEvent.keyDown(field(), { key: "Backspace" });
    fireEvent.keyDown(field(), { key: "Backspace", repeat: true });
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toBeInTheDocument();
    expect(onClose).not.toHaveBeenCalled();
  });

  it("ignores Enter while an input method is composing", () => {
    open();
    fireEvent.keyDown(field(), { key: "Enter", isComposing: true });
    expect(ran).toEqual([]);
  });

  it("goes back to the parent, with what was typed there, when the open page disappears", () => {
    const { onClose, rerender } = open();
    fireEvent.change(field(), { target: { value: "compare" } });
    fireEvent.keyDown(field(), { key: "Enter" });
    const { pages } = setup();
    const { compare: _gone, ...rest } = pages;
    rerender(rest);
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toHaveValue("compare");
    rerender(pages);
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toBeInTheDocument();
    fireEvent.keyDown(field(), { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("names each group for screen readers and keeps the empty message out of the listbox", () => {
    open();
    expect(screen.getByRole("group", { name: "This sample" })).toBeInTheDocument();
    fireEvent.change(field(), { target: { value: "zzz" } });
    expect(screen.getByText("No match.").closest("[role='listbox']")).toBeNull();
  });
});

describe("rootRows ordering", () => {
  it("puts Go to first when the typed text is exactly a sample's name, and caps search rows", () => {
    const { groups, pages } = setup();
    const ids = (q: string): string[] => rootRows(groups, pages, q).flatMap((r) => (r.kind === "item" ? [r.item.id] : []));
    expect(ids("b2")).toEqual(["compare-B2", "go-b2"]);
    const named = groups.map((g) => (g.searchOnly ? { ...g, items: g.items.map((i) => ({ ...i, name: "B2" })) } : g));
    const first = rootRows(named, pages, "B2");
    expect(first[0]).toEqual({ kind: "head", title: "Samples" });
    const many = groups.map((g) => (g.searchOnly ? { ...g, items: Array.from({ length: 50 }, (_, i) => ({ id: `go-${i}`, label: `Go to S${i}` })) } : g));
    expect(rootRows(many, pages, "go to").filter((r) => r.kind === "item")).toHaveLength(20);
  });
});

describe("stepHighlight", () => {
  it("wraps both ways and stays at 0 on an empty list", () => {
    expect(stepHighlight(2, 3, 1)).toBe(0);
    expect(stepHighlight(0, 3, -1)).toBe(2);
    expect(stepHighlight(0, 0, 1)).toBe(0);
  });
});

describe("cycleOption", () => {
  const opts = ["a", "b", "c"].map((v) => ({ value: v, label: v }));
  it("wraps both ways and starts from an end when the value is not an option", () => {
    expect(cycleOption(opts, "c", 1)).toBe("a");
    expect(cycleOption(opts, "a", -1)).toBe("c");
    expect(cycleOption(opts, "zz", 1)).toBe("a");
    expect(cycleOption(opts, "zz", -1)).toBe("c");
  });
});
