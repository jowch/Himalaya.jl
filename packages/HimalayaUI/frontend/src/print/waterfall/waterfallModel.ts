import type { SeriesMember, Trace, AssignmentState } from "../../api";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

/** One indexed-peak bead position on a member's trace. */
export interface WaterfallAnchor {
  /** Effective-peak id (from peak_ids) — maps 1:1 to PlotPeak.id. */
  id: number;
  q: number;
  /** Measured peak intensity (from effective_peaks) — places the bead on the curve. */
  intensity: number | null;
  /** Member dominant phase (drives bead colour); null = unindexed. */
  phase: string | null;
}

/** One member row, fully resolved for rendering (low→high variable order). */
export interface WaterfallRow {
  key: string;
  label: string;
  trace: Trace;
  phase: string | null;
  state: AssignmentState;
  anchors: WaterfallAnchor[];
  bandHeight: number;
  yOffset: number;
}

/** Dominant phase: first confirmed_phases entry, else confirmed_index.phase, else null. */
export function dominantPhase(member: SeriesMember): string | null {
  const snap = member.snapshot;
  if (snap === null) return null;
  if (snap.assignment_state === "form_factor" || snap.assignment_state === "null") return null;
  const cp = snap.confirmed_phases;
  if (cp && cp.length > 0) return cp[0]!.phase;
  return snap.confirmed_index?.phase ?? null;
}

/** Resolve the assignment state, treating missing snapshot or undefined state as "indexed". */
export function resolveState(member: SeriesMember): AssignmentState {
  const snap = member.snapshot;
  if (snap === null) return "indexed";
  return snap.assignment_state ?? "indexed";
}

/** Map API members + a trace lookup → rows in input (display) order. */
export function toWaterfallRows(
  members: SeriesMember[],
  tracesById: Record<number, Trace>,
): WaterfallRow[] {
  return members.map((member) => {
    const snap = member.snapshot;
    const state = resolveState(member);
    const phase = dominantPhase(member);
    const trace = (member.exposure_id != null && tracesById[member.exposure_id]) || EMPTY_TRACE;

    const effectivePeaks = snap?.effective_peaks ?? [];
    const confirmedIndex = snap?.confirmed_index ?? null;

    const peakById = new Map(effectivePeaks.map((p) => [p.id, p]));
    const anchors: WaterfallAnchor[] =
      state === "indexed" && confirmedIndex !== null
        ? confirmedIndex.peak_ids
            .filter((id) => peakById.has(id))
            .map((id) => {
              const p = peakById.get(id)!;
              return { id, q: p.q, intensity: p.intensity, phase };
            })
        : [];

    return {
      key: String(member.id),
      label: member.label_override ?? (member.exposure_id != null ? `exp ${member.exposure_id}` : ""),
      trace,
      phase,
      state,
      anchors,
      bandHeight: member.band_height,
      yOffset: member.y_offset,
    };
  });
}

/** Positive q-extent across all rows, padded ×0.97 / ×1.03. */
export function waterfallQDomain(rows: WaterfallRow[]): [number, number] {
  let lo = Infinity;
  let hi = 0;
  for (const r of rows) {
    for (const q of r.trace.q) {
      if (q > 0 && Number.isFinite(q)) {
        if (q < lo) lo = q;
        if (q > hi) hi = q;
      }
    }
  }
  if (!Number.isFinite(lo) || hi <= 0) return [0.01, 1]; // no positive data fallback
  return [lo * 0.97, hi * 1.03];
}
