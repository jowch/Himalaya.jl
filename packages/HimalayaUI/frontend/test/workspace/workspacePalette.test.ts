import { describe, expect, it, vi } from "vitest";
import type { IndexEntry } from "../../src/api";
import { workspacePalette, type WorkspacePaletteInput } from "../../src/workspace/palette";
import { PLOT_VIEWS } from "../../src/workspace/plotView";

function candidate(id: number, lattice: number): IndexEntry {
  return { id, exposure_id: 1, phase: "Ia3d", basis: 0, score: 1, r_squared: null, lattice_d: lattice, ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [] };
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
    labels: { open: noop, onSample: [], recolour: noop },
    customIndex: noop,
    formFactor: { on: false, toggle: noop },
    views: { options: PLOT_VIEWS, current: "log", set: noop },
    resetZoom: noop,
    comparisons: { shown: true, toggle: noop },
    rescan: undefined,
    problems: undefined,
    removeExperiment: undefined,
    home: noop,
    theme: { current: "system", choose: noop },
    user: { name: "jc", switch: noop },
    help: noop,
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
    expect(items[10]).toMatchObject({ label: "Ia3d 110.0 Å", hint: "Ticked" });
    items[11]!.run!();
    expect(tick).toHaveBeenCalledWith(list[11]);
    expect(groups[0]!.items.find((i) => i.id === "tick")?.hint).toBe("12 candidates");
  });

  it("offers no editing pages or commands before a user is set", () => {
    const { pages, groups } = workspacePalette(
      input({ canEdit: false, candidates: { list: [candidate(1, 100)], ticked: new Set(), tick: vi.fn() }, labels: { open: vi.fn(), onSample: [{ label: "x", name: "x" }], recolour: vi.fn() } }),
    );
    expect(Object.keys(pages).sort()).toEqual(["compare", "theme", "view"]);
    const ids = groups.flatMap((g) => g.items.map((i) => i.id));
    for (const id of ["label", "bridge", "tick", "form-factor", "colours"]) expect(ids).not.toContain(id);
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
});
