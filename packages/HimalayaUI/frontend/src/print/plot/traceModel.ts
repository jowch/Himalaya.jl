import type { Peak, Trace } from "../../api";
import { traceIntensityAt } from "../../lib/plot/traceIntensity";
import type { TraceModel } from "./TracePlot";
import type { PlotPeak } from "./marks/PlotPeaks";

export function toTraceModel(trace: Trace, peaks: readonly Peak[], phase: string | null): TraceModel {
  return {
    trace,
    phase,
    peaks: peaks.map((p): PlotPeak => {
      const intensity = p.intensity != null ? p.intensity : traceIntensityAt(p.q, trace);
      return {
        id: p.id,
        q: p.q,
        source: p.source,
        excluded: p.excluded,
        ...(intensity != null ? { intensity } : {}),
      };
    }),
  };
}
