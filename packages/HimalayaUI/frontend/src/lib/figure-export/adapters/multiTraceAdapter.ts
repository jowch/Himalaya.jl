// multiTraceAdapter.ts — MultiTracePlot state → ExportSpec.
import type { SeriesMember, Trace } from "../../../api";
import type { ExportSpec, LegendRow } from "../types";
import {
  COMPARE_DIMS, COMPARE_PALETTE_LIGHT, ORPHAN_FALLBACK_LIGHT,
  EXPORT_MARGIN,
} from "../presets";
import { buildMultiTraceExportMarks } from "../marks/multiTraceExportMarks";
import { colorFor, ORPHAN_FALLBACK, type GroupingMode } from "../../comparison/coloring";
import { phaseColor } from "../../../phases";

export interface MultiTraceAdapterArgs {
  members: SeriesMember[];          // sorted by display_order
  traces: Map<number, Trace>;
  comparisonTitle: string;
  experimentName?: string;          // omit in /compare/all global scope
  xDomain: [number, number] | null;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  groupingMode: GroupingMode;
  sampleIdFor: (m: SeriesMember) => number | null;
  /** Pre-resolved labels via lib/comparison/labels.resolveDisplayLabels.
   *  Falls back to "Exposure #${exposure_id}" when missing. */
  displayLabelByMemberId?: Map<number, string>;
}

export function buildMultiTraceExportSpec(args: MultiTraceAdapterArgs): ExportSpec {
  const {
    members, traces, comparisonTitle, experimentName,
    xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
  } = args;

  // Critical ordering (spec §MultiTracePlot export "Critical"): compute
  // colors on the UNFILTERED members so defaultDistinct findIndex is stable,
  // then filter null-exposure_id rows from the mark-build pass.
  //
  // colorFor returns the dark-bg ORPHAN_FALLBACK for null-exposure / missing-
  // sample / missing-phase members regardless of which palette we pass — the
  // orphan constant is hardcoded inside coloring.ts. For the export's white
  // background we want the light-bg variant. Swap dark→light here so the
  // marks render orphans with correct contrast AND the legend's color-equality
  // check below detects them.
  const colorByMember = new Map<number, string>();
  for (const m of members) {
    const c = colorFor(m, groupingMode, COMPARE_PALETTE_LIGHT, {
      allMembers: members,
      sampleIdFor,
    });
    colorByMember.set(m.id, c === ORPHAN_FALLBACK ? ORPHAN_FALLBACK_LIGHT : c);
  }

  const filtered = members.filter((m) => m.exposure_id !== null);
  if (filtered.length === 0) {
    throw new Error(
      "buildMultiTraceExportSpec: every member has exposure_id === null; "
      + "export should have been gated by FigureExportControls",
    );
  }

  // Build a default display-label map if the parent didn't provide one. The
  // ComparePage call site DOES provide one; the adapter test calls without.
  const labels = displayLabelByMemberId ?? new Map<number, string>(
    filtered.map((m) => [
      m.id,
      m.label_override ?? `Exposure #${m.exposure_id}`,
    ]),
  );

  // panelHeight for the marks builder = export height minus chrome.
  const panelHeight = COMPARE_DIMS.height - EXPORT_MARGIN.top - EXPORT_MARGIN.bottom;
  const marks = buildMultiTraceExportMarks({
    members: filtered,
    traces,
    displayLabelByMemberId: labels,
    colorByMember,
    showPeakTicks, showPeakLabels,
    panelHeight,
  });

  // Legend per grouping mode.
  const legendRows = buildLegendRows(members, filtered, groupingMode, sampleIdFor, colorByMember);

  const xConfig: Record<string, unknown> = {
    type: "log",
    label: "q",
  };
  if (xDomain) xConfig.domain = xDomain;

  const yConfig: Record<string, unknown> = {
    type: "linear",
    label: null,
    domain: [0, panelHeight], // synthetic — marks pre-place y in pixel space
    axis: null,
  };

  const title: ExportSpec["title"] = experimentName !== undefined
    ? { primary: comparisonTitle, secondary: experimentName }
    : { primary: comparisonTitle };

  return {
    title,
    width: COMPARE_DIMS.width,
    height: COMPARE_DIMS.height,
    plot: {
      marks,
      x: xConfig,
      y: yConfig,
      width: COMPARE_DIMS.width  - EXPORT_MARGIN.left - EXPORT_MARGIN.right,
      height: panelHeight,
      marginLeft: 0,
      marginRight: 0,
      marginTop: 0,
      marginBottom: 0,
    },
    legend: { rows: legendRows },
  };
}

function buildLegendRows(
  unfilteredMembers: SeriesMember[],
  filteredMembers: SeriesMember[],
  mode: GroupingMode,
  sampleIdFor: (m: SeriesMember) => number | null,
  colorByMember: Map<number, string>,
): LegendRow[] {
  if (mode === "distinct") return [];

  const rows: LegendRow[] = [];
  const seen = new Set<string>();
  let orphanPresent = false;

  for (const m of filteredMembers) {
    const color = colorByMember.get(m.id) ?? ORPHAN_FALLBACK_LIGHT;
    if (color === ORPHAN_FALLBACK_LIGHT) {
      orphanPresent = true;
      continue;
    }
    if (mode === "bySample") {
      const sid = sampleIdFor(m);
      const key = `bySample:${sid ?? "null"}`;
      if (seen.has(key)) continue;
      seen.add(key);
      rows.push({
        color,
        symbol: "swatch",
        label: sid !== null ? `Sample ${sid}` : "(unknown sample)",
      });
    } else if (mode === "byPhase") {
      const phase = m.snapshot?.confirmed_index?.phase;
      if (!phase) {
        orphanPresent = true;
        continue;
      }
      if (seen.has(phase)) continue;
      seen.add(phase);
      rows.push({
        color: phaseColor(phase),
        symbol: "swatch",
        label: phase,
      });
    }
  }

  // Note: unfilteredMembers may include null-exposure rows that bind via the
  // bySample sampleIdFor path; those still produce orphan colour and should
  // contribute to the orphan row presence.
  for (const m of unfilteredMembers) {
    if (m.exposure_id === null) {
      orphanPresent = true;
      break;
    }
  }

  if (orphanPresent) {
    rows.push({
      color: ORPHAN_FALLBACK_LIGHT,
      symbol: "swatch",
      label: "unphased / unbound",
    });
  }

  return rows;
}
