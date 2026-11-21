import { describe, expect, it, vi } from "vitest";
import type { IndexEntry } from "../../src/api";
import { workspacePalette, type PaletteLabel, type WorkspacePaletteInput } from "../../src/workspace/palette";
import { PLOT_VIEWS } from "../../src/workspace/plotView";

function candidate(id: number, lattice: number): IndexEntry {
  return { id, exposure_id: 1, phase: "Ia3d", basis: 0, score: 1, r_squared: null, lattice_d: lattice, ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [] };
}

function known(label: string, over: Partial<PaletteLabel> = {}): PaletteLabel {
  return { label, name: label, color: "violet", uses: 1, onSample: undefined, ...over };
}

function input(over: Partial<WorkspacePaletteInput> = {}): WorkspacePaletteInput {
  const noop = vi.fn();
  return {
    modKey: "Ctrl",
    sampleName: "A1",
    canEdit: true,
    nextUnreviewed: noop,
    compare: { samples: [], notes: new Map(), pick: noop },
    bridge: { samples: [], pick: noop, bridged: false, leave: noop },
    candidates: { list: [], ticked: new Set(), tick: noop },
    labels: { known: [], colorFor: () => "plum", toggle: noop, create: noop, recolour: noop },
    customIndex: noop,
    formFactor: { on: false, toggle: noop },
    views: { options: PLOT_VIEWS, current: "log", set: noop },
    resetZoom: noop,
    comparisons: { shown: true, toggle: noop },
    rescan: undefined,
    problems: undefined,
    removeExperiment: undefined,
    settings: undefined,
    home: noop,
    theme: { current: "system", choose: noop },
    user: { name: "jc", switch: noop },
    help: noop,
    history: { undo: undefined, redo: undefined, run: noop },
    samples: [],
    goSample: noop,
    ...over,
  };
}

describe("workspacePalette", () => {
  it("lists every candidate on the tick page, past the nine with number keys, and marks the ticked ones", () => {
    const tick = vi.fn();
    const list = Array.from({ length: 12 }, (_, i) => candidate(i + 1, 100 + i));
    const { pages, groups } = workspacePalette(input({ candidates: { list, ticked: new Set([11]), tick } }));
    const items = pages.tick!.items;
    expect(items).toHaveLength(12);
    expect(items.map((i) => i.shortcut)).toEqual(["1", "2", "3", "4", "5", "6", "7", "8", "9", undefined, undefined, undefined]);
    expect(items[10]).toMatchObject({ label: "Ia3d, a 110.0 Å", hint: "Ticked" });
    items[11]!.run!();
    expect(tick).toHaveBeenCalledWith(list[11]);
    expect(groups[0]!.items.find((i) => i.id === "tick")?.hint).toBe("12 candidates");
  });

  it("offers no editing pages or commands before a user is set", () => {
    const { pages, groups } = workspacePalette(
      input({ canEdit: false, candidates: { list: [candidate(1, 100)], ticked: new Set(), tick: vi.fn() }, labels: { known: [known("x")], colorFor: () => "plum", toggle: vi.fn(), create: vi.fn(), recolour: vi.fn() } }),
    );
    expect(Object.keys(pages).sort()).toEqual(["compare", "theme", "view"]);
    expect(groups.flatMap((g) => g.items.map((i) => i.id))).not.toContain("colours");
    const ids = groups.flatMap((g) => g.items.map((i) => i.id));
    for (const id of ["label", "bridge", "tick", "form-factor"]) expect(ids).not.toContain(id);
  });

  it("names the current view and theme on their rows and pages", () => {
    const { pages, groups } = workspacePalette(input({ views: { options: PLOT_VIEWS, current: "kratky", set: vi.fn() }, theme: { current: "dark", choose: vi.fn() } }));
    const all = groups.flatMap((g) => g.items);
    expect(all.find((i) => i.id === "view")?.hint).toBe("q²I (Kratky)");
    expect(all.find((i) => i.id === "theme")?.hint).toBe("Dark");
    expect(pages.view!.items.filter((i) => i.hint === "Current").map((i) => i.label)).toEqual(["q²I (Kratky)"]);
    expect(pages.theme!.items.filter((i) => i.hint === "Current").map((i) => i.label)).toEqual(["Dark"]);
  });

  it("keeps removing the experiment type-only and marked as danger", () => {
    const { groups } = workspacePalette(input({ removeExperiment: vi.fn() }));
    expect(groups.find((g) => g.title === "Experiment")!.items.find((i) => i.id === "remove-experiment")).toMatchObject({ typeOnly: true, danger: true });
  });
  it("keeps the root order and lists undo and redo only when typed, at the end of App", () => {
    const run = vi.fn();
    const { groups } = workspacePalette(input({ history: { undo: "dropped frame 2 of A1", redo: "ticked Pn3m on A1", run } }));
    expect(groups.map((g) => g.title)).toEqual(workspacePalette(input()).groups.map((g) => g.title));
    expect(groups[0]!.items[0]!.id).toBe("next-unreviewed");
    const app = groups.find((g) => g.title === "App")!.items;
    const undo = app.find((i) => i.id === "undo")!;
    expect(undo).toMatchObject({ label: "Undo: dropped frame 2 of A1", typeOnly: true });
    expect(app.find((i) => i.id === "redo")).toMatchObject({ label: "Redo: ticked Pn3m on A1", typeOnly: true });
    undo.run?.();
    expect(run).toHaveBeenCalledWith("undo");
  });

  it("offers Create for a typed name the experiment does not have yet, which opens a name and colour form", () => {
    const create = vi.fn();
    const { pages } = workspacePalette(input({ labels: { known: [known("retake"), known("Buffer")], colorFor: () => "teal", toggle: vi.fn(), create, recolour: vi.fn() } }));
    expect(pages.label!.extra!("RETAKE")).toEqual([]);
    expect(pages.label!.extra!("")).toEqual([]);
    const [row] = pages.label!.extra!("dehydrated");
    expect(row!.label).toBe("Create “dehydrated”");
    const form = row!.form!();
    expect(form.fields.map((f) => [f.label, f.value])).toEqual([["Name", "dehydrated"], ["Colour", "teal"]]);
    expect(form.fields[1]!.options).toHaveLength(8);
    const name = form.fields[0]!;
    expect(name.check!("")).toBe("Type a name.");
    expect(name.check!("Retake")).toBe("That label exists. Pick it from the list.");
    expect(name.check!("dehydrated")).toBeUndefined();
    form.submit({ name: "  dehydrated ", color: "ochre" });
    expect(create).toHaveBeenCalledWith("dehydrated", "ochre");
  });

  it("labels rows say whether the label is on this sample, and Enter toggles it", () => {
    const toggle = vi.fn();
    const on = known("37C", { onSample: "confirmed" });
    const { pages } = workspacePalette(input({ labels: { known: [on, known("Buffer", { onSample: "suggested" }), known("x")], colorFor: () => "plum", toggle, create: vi.fn(), recolour: vi.fn() } }));
    expect(pages.label!.items.map((i) => i.hint)).toEqual(["On this sample, Enter removes it", "Suggested, Enter confirms it", undefined]);
    pages.label!.items[0]!.run!();
    expect(toggle).toHaveBeenCalledWith(on);
  });

  it("recolours a label through a one-field form that starts on its colour and writes only a change", () => {
    const recolour = vi.fn();
    const { pages } = workspacePalette(input({ labels: { known: [known("37C", { color: "rose", uses: 3 })], colorFor: () => "plum", toggle: vi.fn(), create: vi.fn(), recolour } }));
    const row = pages.colours!.items[0]!;
    expect(row.hint).toBe("3 samples");
    const form = row.form!();
    expect(form.fields.map((f) => f.value)).toEqual(["rose"]);
    form.submit({ color: "rose" });
    expect(recolour).not.toHaveBeenCalled();
    form.submit({ color: "green" });
    expect(recolour).toHaveBeenCalledWith("37C", "green");
  });

  it("lists the experiment name and description on the Settings page only when the user can edit", () => {
    const save = vi.fn();
    expect(workspacePalette(input()).pages.settings).toBeUndefined();
    const { groups, pages } = workspacePalette(input({ settings: { name: "Beamtime 11", description: "", save } }));
    expect(groups.find((g) => g.title === "Experiment")!.items.find((i) => i.id === "settings")?.page).toBe("settings");
    expect(pages.settings!.items.map((i) => [i.label, i.value])).toEqual([
      ["Name", "Beamtime 11"],
      ["Description", ""],
    ]);
  });

  it("saves a trimmed name only when it changed and refuses an empty one", () => {
    const save = vi.fn();
    const { pages } = workspacePalette(input({ settings: { name: "Beamtime 11", description: "", save } }));
    const form = pages.settings!.items[0]!.form!();
    expect(form.fields[0]!.check!("  ")).toBe("Type a name.");
    form.submit({ v: " Beamtime 11 " });
    expect(save).not.toHaveBeenCalled();
    form.submit({ v: " Beamtime 12 " });
    expect(save).toHaveBeenCalledWith({ name: "Beamtime 12" });
  });

  it("trims stray whitespace off a stored name when it is saved unchanged", () => {
    const save = vi.fn();
    const { pages } = workspacePalette(input({ settings: { name: "Beamtime 11 ", description: "", save } }));
    pages.settings!.items[0]!.form!().submit({ v: "Beamtime 11 " });
    expect(save).toHaveBeenCalledWith({ name: "Beamtime 11" });
  });

  it("clears the description to null when it is emptied", () => {
    const save = vi.fn();
    const { pages } = workspacePalette(input({ settings: { name: "B", description: "Lipid screen", save } }));
    pages.settings!.items[1]!.form!().submit({ v: "  " });
    expect(save).toHaveBeenCalledWith({ description: null });
  });
});
