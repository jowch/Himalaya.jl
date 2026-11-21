import { LABEL_COLORS, type CorpusSample, type IndexEntry, type LabelColor } from "../api";
import { cleanLabel, sameLabel } from "../lib/labels";
import type { PlotView, PlotViewOption } from "./plotView";
import { latticeParam } from "./indexing";
import type { WorkspaceTheme } from "./theme";
import { LABEL_SWATCH } from "./ui/labelSwatches";
import type { FieldOption, PaletteForm, PaletteGroup, PaletteItem, PalettePage } from "./ui/palette/model";

export interface PaletteLabel {
  label: string;
  name: string;
  color: LabelColor;
  uses: number;
  onSample: "confirmed" | "suggested" | undefined;
}

export interface WorkspacePaletteInput {
  modKey: string;
  canEdit: boolean;
  nextUnreviewed: () => void;
  compare: { samples: readonly CorpusSample[]; notes: ReadonlyMap<number, string>; pick: (id: number) => void };
  bridge: { samples: readonly CorpusSample[]; pick: (id: number) => void; bridged: boolean; leave: () => void };
  candidates: { list: readonly IndexEntry[]; ticked: ReadonlySet<number>; tick: (ix: IndexEntry) => void };
  labels: {
    known: readonly PaletteLabel[];
    nextColor: LabelColor;
    toggle: (l: PaletteLabel) => void;
    create: (name: string, color: LabelColor) => void;
    recolour: (label: string, color: LabelColor) => void;
  };
  customIndex: () => void;
  formFactor: { on: boolean; toggle: () => void };
  views: { options: readonly PlotViewOption[]; current: PlotView; set: (v: PlotView) => void };
  resetZoom: () => void;
  comparisons: { shown: boolean; toggle: () => void };
  settings: { name: string; description: string; save: (patch: { name?: string; description?: string | null }) => void } | undefined;
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

const COLOR_OPTIONS: readonly FieldOption[] = LABEL_COLORS.map((c) => ({ value: c, label: LABEL_SWATCH[c].name, swatch: LABEL_SWATCH[c].solid }));
const asColor = (v: string): LabelColor => LABEL_COLORS.find((c) => c === v) ?? LABEL_COLORS[0];

function newLabelForm(name: string, p: WorkspacePaletteInput): PaletteForm {
  return {
    title: "New label",
    verb: "create",
    fields: [
      {
        id: "name",
        label: "Name",
        kind: "text",
        value: name,
        check: (v) => {
          const clean = cleanLabel(v);
          if (clean === undefined) return v.trim() === "" ? "Type a name." : "A label is one line of at most 64 characters.";
          return p.labels.known.some((l) => sameLabel(l.label, clean)) ? "That label exists. Pick it from the list." : undefined;
        },
      },
      { id: "color", label: "Colour", kind: "choice", value: p.labels.nextColor, options: COLOR_OPTIONS },
    ],
    submit: (v) => p.labels.create(cleanLabel(v.name ?? "")!, asColor(v.color ?? "")),
  };
}

function recolourForm(l: PaletteLabel, p: WorkspacePaletteInput): PaletteForm {
  return {
    title: l.name,
    verb: "apply",
    fields: [{ id: "color", label: "Colour", kind: "choice", value: l.color, options: COLOR_OPTIONS }],
    submit: (v) => {
      const c = asColor(v.color ?? "");
      if (c !== l.color) p.labels.recolour(l.label, c);
    },
  };
}

function textForm(title: string, value: string, check: (v: string) => string | undefined, save: (v: string) => void): PaletteForm {
  return {
    title,
    verb: "save",
    fields: [{ id: "v", label: title, kind: "text", value, check }],
    submit: (v) => {
      const next = (v.v ?? "").trim();
      if (next !== value.trim()) save(next);
    },
  };
}

function settingsPage(s: NonNullable<WorkspacePaletteInput["settings"]>): PalettePage {
  return {
    title: "Settings",
    placeholder: "Find a setting",
    items: [
      {
        id: "setting-name",
        label: "Name",
        value: s.name,
        form: () => textForm("Name", s.name, (v) => (v.trim() ? undefined : "Type a name."), (name) => s.save({ name })),
      },
      {
        id: "setting-description",
        label: "Description",
        value: s.description,
        form: () => textForm("Description", s.description, () => undefined, (d) => s.save({ description: d || null })),
      },
    ],
  };
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
  if (p.settings) pages.settings = settingsPage(p.settings);
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
    pages.label = {
      title: "Label",
      arg: "Name",
      placeholder: "Pick or type a new label",
      deep: true,
      items: p.labels.known.map((l) => ({
        id: `label-${l.label}`,
        label: l.name,
        swatch: LABEL_SWATCH[l.color].solid,
        hint: l.onSample === "confirmed" ? "On this sample, Enter removes it" : l.onSample === "suggested" ? "Suggested, Enter confirms it" : undefined,
        run: () => p.labels.toggle(l),
      })),
      extra: (q) => {
        const name = cleanLabel(q);
        if (name === undefined || p.labels.known.some((l) => sameLabel(l.label, name))) return [];
        return [{ id: "label-create", label: `Create “${name}”`, form: () => newLabelForm(name, p) }];
      },
    };
    pages.colours = {
      title: "Label colour",
      arg: "Label",
      placeholder: "Pick a label",
      items: p.labels.known.map((l) => ({
        id: `colour-${l.label}`,
        label: l.name,
        swatch: LABEL_SWATCH[l.color].solid,
        hint: `${l.uses} sample${l.uses === 1 ? "" : "s"}`,
        form: () => recolourForm(l, p),
      })),
    };
  }

  const n = p.candidates.list.length;
  const sample: PaletteItem[] = [
    { id: "next-unreviewed", label: "Next unreviewed sample", shortcut: "Z", run: p.nextUnreviewed },
    ...(p.canEdit ? [{ id: "label", label: "Label", page: "label" }] : []),
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
    ...(pages.settings ? [{ id: "settings", label: "Settings", hint: "Name, description", page: "settings" }] : []),
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
