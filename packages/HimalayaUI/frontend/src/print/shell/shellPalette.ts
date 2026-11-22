import { appItems, themePage, type AppActions } from "../../workspace/palette";
import type { PaletteGroup, PaletteItem, PalettePage } from "../../workspace/ui/palette/model";

export interface ShellPaletteInput extends AppActions {
  pathname: string;
  experiments: readonly { id: number; name: string }[];
  samples: readonly { id: number; name: string; experiment: string; href: string }[];
  go: (path: string) => void;
}

export function shellPalette(p: ShellPaletteInput): { groups: PaletteGroup[]; pages: Record<string, PalettePage> } {
  const experiments: PaletteItem[] = [
    ...(p.pathname === "/experiments/new" ? [] : [{ id: "new-experiment", label: "New experiment", run: () => p.go("/experiments/new") }]),
    ...(p.pathname === "/experiments" ? [] : [{ id: "home", label: "All experiments", run: () => p.go("/experiments") }]),
    ...p.experiments.map((e) => ({ id: `open-${e.id}`, label: `Open ${e.name}`, name: e.name, run: () => p.go(`/experiments/${e.id}`) })),
  ];
  return {
    groups: [
      { title: "Experiments", items: experiments },
      { title: "App", items: appItems(p) },
      { title: "Samples", searchOnly: true, items: p.samples.map((s) => ({ id: `go-${s.id}`, label: `Go to ${s.name}`, name: s.name, hint: s.experiment, run: () => p.go(s.href) })) },
    ],
    pages: { theme: themePage(p.theme) },
  };
}
