import { appItems, goToSampleItem, openExperimentItem, themePage, type AppActions } from "../workspace/palette";
import type { CorpusExperiment, CorpusSampleRow } from "../hooks/useCorpusPalette";
import type { PaletteGroup, PaletteItem, PalettePage } from "../ui/palette/model";

export interface ShellPaletteInput extends AppActions {
  pathname: string;
  experiments: readonly CorpusExperiment[];
  samples: readonly CorpusSampleRow[];
  go: (path: string) => void;
}

export function shellPalette(p: ShellPaletteInput): { groups: PaletteGroup[]; pages: Record<string, PalettePage> } {
  const experiments: PaletteItem[] = [
    ...(p.pathname === "/experiments/new" ? [] : [{ id: "new-experiment", label: "New experiment", run: () => p.go("/experiments/new") }]),
    ...(p.pathname === "/experiments" ? [] : [{ id: "home", label: "All experiments", run: () => p.go("/experiments") }]),
    ...p.experiments.map((e) => openExperimentItem(e, p.go)),
  ];
  return {
    groups: [
      { title: "Experiments", items: experiments },
      { title: "App", items: appItems(p) },
      { title: "Samples", searchOnly: true, items: p.samples.map((s) => goToSampleItem(s, p.go)) },
    ],
    pages: { theme: themePage(p.theme) },
  };
}
