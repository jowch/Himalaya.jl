/**
 * MemberHeatmapLayer — per-member mark factory for the **heatmap**
 * representation of the shared `MultiTracePlot` render core (#208).
 *
 * **Shape:** This is NOT a React component. Same contract as
 * `MemberTraceLayer`: given a comparison member + its live `(q, I)` trace +
 * its y-band, return an array of Observable Plot marks the parent plot
 * composes into a single `<Plot>`. One shared plot, one shared q-scale; brush
 * + zoom apply uniformly across all rows.
 *
 * Visual: each member is one horizontal row of `NB` rectangles spanning the
 * full q-domain. Cell fill is the member's grouping/phase color mixed against
 * the paper plate, with mix percentage driven by the binned intensity (peak
 * features map to opaque cells; the slow form-factor decay stays transparent
 * so the heatmap reads "peaks, not slope"). Matches the mockup's
 * `renderHeatmap` (`series-builder.html:758-823`).
 *
 * Y-positioning uses the same y-band envelope the waterfall uses — the only
 * thing that changes between representations is the per-band mark vocabulary.
 *
 * **Peaks-only intensity.** The waterfall renders the full trace; the heatmap
 * intentionally drops the form-factor decay (peaks are the signal here, not
 * the background slope). For each q-bin we take the maximum intensity *above
 * a local background baseline* in that bin — robust enough for the SAXS
 * peak-on-decay shape without re-fitting peaks. See `binPeaksOnly` below.
 */
import * as Plot from "@observablehq/plot";
import type { Trace, SeriesMember } from "../api";
import { phaseColor } from "../phases";
import {
  colorFor,
  COMPARE_PALETTE,
  type GroupingMode,
} from "../lib/comparison/coloring";

/**
 * Paper plate colour (`--plate`) — concrete value so OKLCH `color-mix` calls
 * have a non-`var()` second argument. Mirrors `series-builder.html` line 516
 * (the mockup's `COL.plate`). Pinned in `DESIGN.md` §2 Neutrals.
 */
const PLATE_OKLCH = "oklch(0.992 0.004 90)";

/**
 * Number of q-bins per heatmap row. 200 matches the mockup's `NB = 200` and
 * lands near-screen-resolution at typical builder widths (~6 px per cell at
 * a 1200 px plot interior). Higher counts add render cost without visible
 * resolution; lower counts pixelate peak edges.
 */
export const HEATMAP_BIN_COUNT = 200;

/**
 * Minimum and maximum percent of the cell fill that is "ink" — the rest is
 * paper plate. `MIN_PCT` keeps even empty bins faintly tinted so each row
 * reads as a continuous strip (the mockup uses `6%`); `MAX_PCT` ceilings the
 * saturation so adjacent peaks don't crush into solid blocks (`100%`).
 */
const MIN_PCT = 6;
const MAX_PCT = 100;

/**
 * Power-law contrast curve. `< 1` boosts mid-intensity bins (peak shoulders
 * become visible); `1` is linear. 0.62 matches the mockup. Keep in lockstep
 * with the mockup so the live render and the design reference agree on
 * "what a peak looks like".
 */
const CONTRAST_GAMMA = 0.62;

export interface MemberHeatmapMarksProps {
  member: SeriesMember;
  trace: Trace | undefined;
  /** Numeric y-band envelope `[topPx, bottomPx]` for this row. */
  yBand: [number, number];
  /** Visible q-domain `[qMin, qMax]` (full-range when `null` upstream). */
  qDomain: [number, number];
  /**
   * Grouping context (Phase 9 contract, same as `MemberTraceLayer`). Drives
   * the per-row hue:
   *   - `byPhase` + indexed → phase hue from `phases.ts`
   *   - `bySample`/`distinct` → entry from `COMPARE_PALETTE` via `colorFor`
   *   - omitted → fallback ink (legacy callers without grouping context)
   *
   * `member.color_override` always wins (the resolution order lives inside
   * `colorFor`). When `groupingMode` is omitted, the layer falls back to
   * `phaseColor` for indexed members and a neutral ink for the rest.
   */
  groupingMode?: GroupingMode;
  allMembers?: ReadonlyArray<SeriesMember>;
  sampleIdFor?: (m: SeriesMember) => number | null;
}

/**
 * Build the Observable Plot marks for one heatmap row. Pure function; safe
 * to call inside a render. Emits one `Plot.rect` per cell so the cells share
 * the plot's q-scale (clipping + brush + zoom apply uniformly across the
 * whole stack of rows, exactly like the waterfall traces).
 */
export function buildMemberHeatmapMarks(props: MemberHeatmapMarksProps): unknown[] {
  const marks: unknown[] = [];
  const { member, trace, yBand, qDomain } = props;
  if (!trace || trace.q.length === 0) return marks;

  const [yTop, yBot] = yBand;
  // Inset the row a couple of px so adjacent rows don't visually fuse.
  const rowInset = 2;
  const y1 = yTop + rowInset;
  const y2 = Math.max(y1 + 1, yBot - rowInset);

  // Per-member fill colour — same lookup the waterfall uses.
  const fillBase = resolveHeatmapFill(props);

  // Bin the trace's peak signal across the visible q-domain.
  const bins = binPeaksOnly(trace.q, trace.I, qDomain[0], qDomain[1], HEATMAP_BIN_COUNT);

  // Single-row max for the per-row normaliser (so a quiet sample's strongest
  // peak still reaches MAX_PCT). A cross-row normaliser is plausible but
  // makes faint samples disappear; per-row matches the mockup's `gmax` per
  // sample → per-row band ratio behaviour at the visual level.
  let rowMax = 0;
  for (const v of bins) if (v > rowMax) rowMax = v;
  const norm = rowMax > 0 ? rowMax : 1;

  const cells = bins.map((v, i) => {
    const f0 = i / HEATMAP_BIN_COUNT;
    const f1 = (i + 1) / HEATMAP_BIN_COUNT;
    const q0 = qFromFrac(f0, qDomain[0], qDomain[1]);
    const q1 = qFromFrac(f1, qDomain[0], qDomain[1]);
    const frac = Math.pow(Math.max(0, v) / norm, CONTRAST_GAMMA);
    const pct = MIN_PCT + frac * (MAX_PCT - MIN_PCT);
    return {
      x1: q0,
      x2: q1,
      // The mockup mixes in OKLab; we use OKLCH-compatible mix which renders
      // in the same colour space for `oklch(...)` inputs.
      fill: `color-mix(in oklab, ${fillBase} ${pct.toFixed(0)}%, ${PLATE_OKLCH})`,
      memberId: member.id,
    };
  });

  marks.push(
    Plot.rect(cells, {
      x1: "x1",
      x2: "x2",
      y1: y1,
      y2: y2,
      fill: (d: unknown) => (d as { fill: string }).fill,
      // No hairline — adjacent cells should read as a continuous gradient.
      stroke: null,
    }),
  );

  return marks;
}

/**
 * Resolve the per-row "ink" colour. Mirrors `MemberTraceLayer`'s grouping
 * contract: a member always renders in its chosen group identity (so byPhase
 * vs bySample swap recolors the heatmap exactly the same way it recolors
 * the waterfall).
 */
function resolveHeatmapFill(props: MemberHeatmapMarksProps): string {
  const { member, groupingMode, allMembers, sampleIdFor } = props;
  if (groupingMode !== undefined && allMembers !== undefined && sampleIdFor !== undefined) {
    return colorFor(member, groupingMode, COMPARE_PALETTE, {
      allMembers,
      sampleIdFor,
    });
  }
  // Legacy fallback — phase hue if indexed, neutral ink otherwise.
  const phase = member.snapshot?.confirmed_index?.phase;
  if (phase) return phaseColor(phase);
  return "var(--color-fg)";
}

/**
 * Bin the trace into `n` equal q-bins across `[qMin, qMax]`. Each cell value
 * is the max intensity in that bin minus the row's running baseline (the
 * "peaks, not slope" requirement: a smooth decay should fade to background,
 * a sharp peak should pop). Out-of-domain bins are zero.
 *
 * Implementation: walk the trace once, accumulate per-bin max, then subtract
 * a low-percentile (here min) per bin. For SAXS-shaped traces the smooth
 * decay's `min` is close to the true baseline, so peak amplitude survives
 * subtraction. This is intentionally coarser than the mockup's analytical
 * `peaksOnly` (which knows the form factor) because we don't have a
 * peak-fit model here.
 */
export function binPeaksOnly(
  qs: ReadonlyArray<number>,
  Is: ReadonlyArray<number>,
  qMin: number,
  qMax: number,
  n: number,
): number[] {
  if (qs.length === 0 || qMax <= qMin || n <= 0) return new Array(n).fill(0);
  const out = new Array<number>(n).fill(0);
  // Walk and keep per-bin max.
  for (let i = 0; i < qs.length; i++) {
    const q = qs[i]!;
    if (q < qMin || q > qMax) continue;
    const frac = (q - qMin) / (qMax - qMin);
    const bin = Math.min(n - 1, Math.max(0, Math.floor(frac * n)));
    const v = Is[i] ?? 0;
    if (v > out[bin]!) out[bin] = v;
  }
  // Subtract a smooth running baseline (low-percentile pooled across a sliding
  // window). For SAXS, a rolling-min approximates the decay envelope cheaply.
  const W = Math.max(3, Math.floor(n / 12)); // ~16-bin window at n=200
  const baseline = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) {
    let m = Infinity;
    const lo = Math.max(0, i - W);
    const hi = Math.min(n - 1, i + W);
    for (let j = lo; j <= hi; j++) {
      const v = out[j]!;
      if (v > 0 && v < m) m = v;
    }
    baseline[i] = Number.isFinite(m) ? m : 0;
  }
  for (let i = 0; i < n; i++) {
    const v = out[i]! - baseline[i]!;
    out[i] = v > 0 ? v : 0;
  }
  return out;
}

/** Map a normalised fraction back to q (linear-in-fraction over the domain). */
function qFromFrac(f: number, qMin: number, qMax: number): number {
  return qMin + f * (qMax - qMin);
}
