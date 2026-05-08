// traceAdapter.ts — TraceViewer state → ExportSpec.
import type { Trace, Peak, IndexEntry } from "../../../api";
import type { ExportSpec, LegendRow } from "../types";
import { TRACE_DIMS, LIGHT_PALETTE, EXPORT_MARGIN } from "../presets";
import { buildTraceExportMarks } from "../marks/traceExportMarks";
import { phaseColor } from "../../../phases";

export interface TraceAdapterArgs {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
  experimentName: string;
  sampleName: string;
  exposureLabel: string;
  xDomain: [number, number] | null;
  yDomain: [number, number] | null;
  xType: "log" | "linear";
  /** Display string for the q axis label. Falls back to "Å⁻¹". */
  qUnits?: string;
}

export function buildTraceExportSpec(args: TraceAdapterArgs): ExportSpec {
  const {
    trace, peaks, activeGroupIndices,
    experimentName, sampleName, exposureLabel,
    xDomain, yDomain, xType, qUnits,
  } = args;

  const marks = buildTraceExportMarks({ trace, peaks, activeGroupIndices });

  const xConfig: Record<string, unknown> = {
    type: xType,
    label: `q (${qUnits ?? "Å⁻¹"})`,
  };
  if (xDomain) xConfig.domain = xDomain;

  const yConfig: Record<string, unknown> = {
    type: "log",
    label: "I",
  };
  if (yDomain) yConfig.domain = yDomain;

  const plotW = TRACE_DIMS.width  - EXPORT_MARGIN.left - EXPORT_MARGIN.right;
  const plotH = TRACE_DIMS.height - EXPORT_MARGIN.top  - EXPORT_MARGIN.bottom;

  const legend: LegendRow[] = [];

  // Peak source rows.
  legend.push({ color: LIGHT_PALETTE.peakAuto, symbol: "swatch", label: "auto peak" });
  if (peaks.some((p) => p.source === "manual")) {
    legend.push({ color: LIGHT_PALETTE.peakManual, symbol: "swatch", label: "manual peak" });
  }
  if (peaks.some((p) => p.excluded)) {
    legend.push({ color: LIGHT_PALETTE.peakExcluded, symbol: "swatch", label: "excluded auto peak" });
  }

  // Phase rows — one per active-group index phase.
  const seenPhases = new Set<string>();
  for (const idx of activeGroupIndices) {
    if (seenPhases.has(idx.phase)) continue;
    seenPhases.add(idx.phase);
    legend.push({
      color: phaseColor(idx.phase),
      symbol: "line",
      label: `predicted ${idx.phase}`,
    });
  }

  return {
    title: { primary: `${experimentName} · ${sampleName} · ${exposureLabel}` },
    width: TRACE_DIMS.width,
    height: TRACE_DIMS.height,
    plot: {
      marks,
      x: xConfig,
      y: yConfig,
      width: plotW,
      height: plotH,
      marginLeft: 0,
      marginRight: 0,
      marginTop: 0,
      marginBottom: 0,
    },
    legend: { rows: legend },
  };
}
