import type { CorpusSample, IndexEntry } from "../api";
import type { PlotView, PlotViewOption } from "./plotView";
import { latticeParam } from "./indexing";
import type { WorkspaceTheme } from "./theme";
import type { PaletteGroup, PaletteItem, PalettePage } from "./ui/palette/model";

export interface WorkspacePaletteInput {
  modKey: string;
  canEdit: boolean;
  nextUnreviewed: () => void;
  compare: { samples: readonly CorpusSample[]; notes: ReadonlyMap<number, string>; pick: (id: number) => void };
  bridge: { samples: readonly CorpusSample[]; pick: (id: number) => void; bridged: boolean; leave: () => void };
  candidates: { list: readonly IndexEntry[]; ticked: ReadonlySet<number>; tick: (ix: IndexEntry) => void };
  labels: { open: () => void; onSample: readonly { label: string; name: string }[]; recolour: (label: string) => void };
  customIndex: () => void;
  formFactor: { on: boolean; toggle: () => void };
  views: { options: readonly PlotViewOption[]; current: PlotView; set: (v: PlotView) => void };
  resetZoom: () => void;
  comparisons: { shown: boolean; toggle: () => void };
  rescan: (() => void) | undefined;
  problems: (() => void) | undefined;
  removeExperiment: (() => void) | undefined;
  home: () => void;
  theme: { current: WorkspaceTheme; choose: (t: WorkspaceTheme) => void };
  user: { name: string | undefined; switch: () => void };
  help: () => void;
  samples: readonly CorpusSample[];
  goSample: (id: number) => void;
}

const THEMES: readonly { theme: WorkspaceTheme; name: string }[] = [
  { theme: "light", name: "Light" },
  { theme: "dark", name: "Dark" },
  { theme: "system", name: "Follow the system" },
];

export function candidateName(ix: Pick<IndexEntry, "phase" | "lattice_d">): string {
  return ix.lattice_d !== null ? `${ix.phase}, ${latticeParam(ix.phase)} ${ix.lattice_d.toFixed(1)} Å` : ix.phase;
}

export function workspacePalette(p: WorkspacePaletteInput): { groups: PaletteGroup[]; pages: Record<string, PalettePage> } {
  const pages: Record<string, PalettePage> = {
    compare: {
      title: "Compare with",
      arg: "Sample",
      deep: true,
      empty: "No other samples.",
      items: p.compare.samples.map((s) => ({ id: `compare-${s.id}`, label: s.name, hint: p.compare.notes.get(s.id), run: () => p.compare.pick(s.id) })),
    },
    view: {
      title: "Plot view",
      arg: "View",
      deep: true,
      items: p.views.options.map((o) => ({ id: `view-${o.view}`, label: o.name, hint: o.view === p.views.current ? "Current" : undefined, run: () => p.views.set(o.view) })),
    },
    theme: {
      title: "Theme",
      arg: "Theme",
      items: THEMES.map((t) => ({ id: `theme-${t.theme}`, label: t.name, hint: t.theme === p.theme.current ? "Current" : undefined, run: () => p.theme.choose(t.theme) })),
    },
  };
  if (p.canEdit) {
    pages.bridge = {
      title: "Bridge to",
      arg: "Sample",
      deep: true,
      empty: "No samples left to bridge to.",
      items: p.bridge.samples.map((s) => ({ id: `bridge-${s.id}`, label: s.name, run: () => p.bridge.pick(s.id) })),
    };
    if (p.candidates.list.length > 0) {
      pages.tick = {
        title: "Tick",
        arg: "Candidate",
        placeholder: "Phase or lattice, such as Ia3d 171",
        deep: true,
        items: p.candidates.list.map((ix, i) => ({
          id: `tick-${ix.id}`,
          label: candidateName(ix),
          hint: p.candidates.ticked.has(ix.id) ? "Ticked" : undefined,
          shortcut: i < 9 ? String(i + 1) : undefined,
          run: () => p.candidates.tick(ix),
        })),
      };
    }
    if (p.labels.onSample.length > 0) {
      pages.colours = {
        title: "Label colour",
        arg: "Label",
        placeholder: "Pick a label",
        items: p.labels.onSample.map((l) => ({ id: `colour-${l.label}`, label: l.name, run: () => p.labels.recolour(l.label) })),
      };
    }
  }

  const n = p.candidates.list.length;
  const sample: PaletteItem[] = [
    { id: "next-unreviewed", label: "Next unreviewed sample", shortcut: "Z", run: p.nextUnreviewed },
    ...(p.canEdit ? [{ id: "label", label: "Label", shortcut: "L", run: p.labels.open }] : []),
    { id: "compare", label: "Compare with a sample", shortcut: "C", page: "compare" },
    ...(p.canEdit ? [{ id: "bridge", label: "Bridge to a sample", page: "bridge" }] : []),
    ...(p.canEdit && p.bridge.bridged ? [{ id: "unbridge", label: "Take this sample out of its bridges", run: p.bridge.leave }] : []),
    ...(pages.tick ? [{ id: "tick", label: "Tick a candidate", hint: `${n} candidate${n === 1 ? "" : "s"}`, page: "tick" }] : []),
    { id: "custom", label: "Custom index", shortcut: "I", run: p.customIndex },
    ...(p.canEdit ? [{ id: "form-factor", label: p.formFactor.on ? "Clear no phase (form factor)" : "Mark no phase (form factor)", shortcut: "0", run: p.formFactor.toggle }] : []),
  ];
  const current = p.views.options.find((o) => o.view === p.views.current)?.name;
  const plot: PaletteItem[] = [
    { id: "view", label: "Plot view", hint: current, page: "view" },
    { id: "zoom", label: "Reset zoom", shortcut: `${p.modKey} 0`, run: p.resetZoom },
    { id: "refs", label: p.comparisons.shown ? "Hide comparisons" : "Show comparisons", run: p.comparisons.toggle },
  ];
  const experiment: PaletteItem[] = [
    ...(p.rescan ? [{ id: "rescan", label: "Rescan the data folder", run: p.rescan }] : []),
    ...(p.problems ? [{ id: "problems", label: "Show ingest problems", run: p.problems }] : []),
    { id: "home", label: "Back to all experiments", run: p.home },
    ...(p.removeExperiment ? [{ id: "remove-experiment", label: "Remove this experiment", danger: true, typeOnly: true, run: p.removeExperiment }] : []),
  ];
  const themeName = THEMES.find((t) => t.theme === p.theme.current)?.name;
  const app: PaletteItem[] = [
    ...(pages.colours ? [{ id: "colours", label: "Label colours", page: "colours" }] : []),
    { id: "theme", label: "Theme", hint: themeName, page: "theme" },
    { id: "switch-user", label: p.user.name ? "Switch user" : "Sign in", hint: p.user.name, run: p.user.switch },
    { id: "help", label: "Keyboard shortcuts", shortcut: "?", run: p.help },
  ];
  const groups: PaletteGroup[] = [
    { title: "This sample", items: sample },
    { title: "Plot", items: plot },
    { title: "Experiment", items: experiment },
    { title: "App", items: app },
    { title: "Samples", searchOnly: true, items: p.samples.map((s) => ({ id: `go-${s.id}`, label: `Go to ${s.name}`, name: s.name, run: () => p.goSample(s.id) })) },
  ];
  return { groups, pages };
}
