import { useCallback, useState } from "react";
import type { ScaleType } from "../print/plot/projection";

export type PlotView = "log" | "linear" | "loglog" | "kratky";

export interface PlotViewOption {
  view: PlotView;
  label: string;
  title: string;
}

export const PLOT_VIEWS: readonly PlotViewOption[] = [
  { view: "log", label: "Log", title: "Log intensity against q" },
  { view: "linear", label: "Linear", title: "Linear intensity against q" },
  { view: "loglog", label: "Log-log", title: "Log intensity against log q" },
  { view: "kratky", label: "q²I", title: "q² × intensity against q (Kratky). Lifts weak high-order peaks" },
];

const KEY = "himalaya-ui:workspace-plot-view";

export function viewAxes(view: PlotView): { xType: ScaleType; yType: ScaleType } {
  return {
    xType: view === "loglog" ? "log" : "linear",
    yType: view === "log" || view === "loglog" ? "log" : "linear",
  };
}

// physics: a powder pattern loses intensity with q by a Lorentz factor of 1/q², so
// physics: q²·I (the Kratky axis) restores the relative weight of higher-order
// physics: reflections and makes weak √14 or √16 peaks visible. It is a display
// physics: transform only: peak finding always runs on the measured I.
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
