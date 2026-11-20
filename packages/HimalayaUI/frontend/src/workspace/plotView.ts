import { useCallback, useState } from "react";
import type { ScaleType } from "../print/plot/projection";

export type PlotView = "log" | "linear" | "loglog" | "kratky";

export interface PlotViewOption {
  view: PlotView;
  label: string;
  name: string;
  title: string;
}

export const PLOT_VIEWS: readonly PlotViewOption[] = [
  { view: "log", label: "Log", name: "Log", title: "Log intensity against q" },
  { view: "linear", label: "Linear", name: "Linear", title: "Linear intensity against q" },
  { view: "loglog", label: "Log-log", name: "Log-log", title: "Log intensity against log q" },
  { view: "kratky", label: "q²I", name: "q²I (Kratky)", title: "q² × intensity against q (Kratky)" },
];

const KEY = "himalaya-ui:workspace-plot-view";

export function viewAxes(view: PlotView): { xType: ScaleType; yType: ScaleType } {
  return {
    xType: view === "loglog" ? "log" : "linear",
    yType: view === "log" || view === "loglog" ? "log" : "linear",
  };
}

// physics: the integrated intensity of a powder reflection falls with q by a Lorentz
// physics: factor of 1/q², so q²·I (the Kratky axis) restores the relative weight of
// physics: higher-order reflections. A flat solvent background also grows as q², so on
// physics: SSRL traces the high-q end can sit as high as the strongest peak. It is a
// physics: display transform only: peak finding always runs on the measured I.
export function viewIntensity(view: PlotView, q: number, I: number): number {
  return view === "kratky" ? q * q * I : I;
}

export function applyView<T extends { q: number[]; I: number[]; sigma: number[] }>(view: PlotView, trace: T): T {
  if (view !== "kratky") return trace;
  return {
    ...trace,
    I: trace.I.map((v, i) => viewIntensity(view, trace.q[i]!, v)),
    sigma: trace.sigma.map((v, i) => viewIntensity(view, trace.q[i]!, v)),
  };
}

export function viewLabel(view: PlotView): string {
  return view === "kratky" ? "q²·I (a.u.)" : "intensity (a.u.)";
}

export function nextView(view: PlotView, dir: 1 | -1): PlotView {
  const i = PLOT_VIEWS.findIndex((v) => v.view === view);
  const n = PLOT_VIEWS.length;
  return PLOT_VIEWS[(((i + dir) % n) + n) % n]!.view;
}

function read(): PlotView {
  try {
    const v = window.localStorage.getItem(KEY);
    return PLOT_VIEWS.some((o) => o.view === v) ? (v as PlotView) : "log";
  } catch {
    return "log";
  }
}

export function usePlotView(): [PlotView, (next: PlotView) => void] {
  const [view, setView] = useState<PlotView>(read);
  const choose = useCallback((next: PlotView) => {
    setView(next);
    try {
      window.localStorage.setItem(KEY, next);
    } catch {
      return;
    }
  }, []);
  return [view, choose];
}
