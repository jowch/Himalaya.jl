/**
 * 1-D label dodge layout for the Compare page peak labels
 * (Plan §Phase 8, Task 8.2; spec §Label placement).
 *
 * Labels render vertically above triangle peak markers. When peaks are
 * crowded, naive direct-above placement causes label collisions. This
 * function operates in *pixel space* on a list of labeled-peak positions,
 * spreads them horizontally to enforce a minimum pixel gap (= label width),
 * and emits per-label coordinates plus the original triangle position so a
 * leader line can connect text → triangle for displaced labels.
 *
 * Algorithm: a single left-to-right sweep that nudges each label rightward
 * if it would collide with the previous label, then a right-to-left
 * smoothing pass to keep collisions resolved while pulling labels back
 * toward their natural q (so the dodge is symmetric for clusters and stays
 * out of the way for sparse layouts).
 *
 * Why pure pixel-space + a callback pair (toPx / fromPx) rather than
 * data-space arithmetic: the dodge needs comparable horizontal extents.
 * In log-q the gap between q=0.30 and q=0.32 in pixels depends on the
 * scale and the plot's domain — only pixel-space dodging gives consistent
 * label spacing.
 *
 * Output:
 *   - `qPeak`, `yPeak` — the original triangle position (anchor for the
 *     leader line). Equal to the input `q`/`y`.
 *   - `qLabel`, `yLabel` — the rendered label position. `qLabel === qPeak`
 *     when no dodge was needed; differs when the layout pushed it sideways.
 *   - `label` — the rendered text (q to 3 sig digits, matching v6.1 format).
 *   - `peakId` — passed through for downstream callers (e.g. tooltip lookup).
 *
 * No-op for empty input or a single peak (returns it unchanged at qLabel
 * = qPeak).
 */

/** Pixel offset above the peak triangle at which labels render. */
const LABEL_OFFSET_PX = 12;

export interface LabelInputPeak {
  q: number;
  y: number;
  peakId: number;
}

export interface LabelDodgeOpts {
  /** q → pixel x for the active x-scale. */
  toPx: (q: number) => number;
  /** pixel x → q. Used to convert dodged pixel positions back to data space. */
  fromPx: (px: number) => number;
  /**
   * Minimum horizontal extent each label requires (estimate). Pairs of
   * labels closer than this in pixel space get pushed apart.
   */
  labelWidthPx: number;
}

export interface DodgedLabel {
  /** Original triangle q (leader-line anchor). */
  qPeak: number;
  /** Original triangle y in pixel space (leader-line anchor). */
  yPeak: number;
  /** Rendered label q. Equal to qPeak when no dodge was needed. */
  qLabel: number;
  /** Rendered label y in pixel space (above the triangle). */
  yLabel: number;
  /** Display text for the label. */
  label: string;
  /** Pass-through peak identifier. */
  peakId: number;
}

export function layoutPeakLabels(
  peaks: LabelInputPeak[],
  opts: LabelDodgeOpts,
): DodgedLabel[] {
  if (peaks.length === 0) return [];

  const { toPx, fromPx, labelWidthPx } = opts;

  // Sort by natural pixel position so the sweep makes sense; remember the
  // original input index so we can restore order in the output.
  const indexed = peaks.map((p, idx) => ({
    idx,
    p,
    naturalPx: toPx(p.q),
  }));
  indexed.sort((a, b) => a.naturalPx - b.naturalPx);

  // Sweep left → right: each label must be at least `labelWidthPx` to the
  // right of the previous one. This produces a candidate layout that
  // respects the gap but shifts the cluster rightward (asymmetric).
  const naturalPx = indexed.map((it) => it.naturalPx);
  const dodgedPx: number[] = new Array(indexed.length);
  dodgedPx[0] = naturalPx[0]!;
  for (let i = 1; i < indexed.length; i++) {
    const minPx = dodgedPx[i - 1]! + labelWidthPx;
    dodgedPx[i] = Math.max(naturalPx[i]!, minPx);
  }
  // Recenter: shift the whole layout so the average dodged position equals
  // the average natural position. Keeps the cluster symmetric around its
  // original mean (one label nudges left, the next right, etc.).
  let naturalSum = 0;
  let dodgedSum = 0;
  for (let i = 0; i < indexed.length; i++) {
    naturalSum += naturalPx[i]!;
    dodgedSum += dodgedPx[i]!;
  }
  const shift = (naturalSum - dodgedSum) / indexed.length;
  if (shift !== 0) {
    for (let i = 0; i < indexed.length; i++) dodgedPx[i] = dodgedPx[i]! + shift;
  }
  // The recenter may have introduced a violation only if the left→right
  // sweep wasn't already gap-correct; since the sweep guarantees gap
  // satisfaction and a uniform shift preserves gaps, no re-fix needed.
  // For sparse layouts (no collisions), all dodgedPx[i] == naturalPx[i]
  // and the shift is zero — labels stay put.
  // Snap labels back to natural when the shift didn't change anything
  // (within float precision) so qLabel === qPeak holds exactly for the
  // sparse case (relied on by the leader-line filter).
  for (let i = 0; i < indexed.length; i++) {
    if (Math.abs(dodgedPx[i]! - naturalPx[i]!) < 1e-9) {
      dodgedPx[i] = naturalPx[i]!;
    }
  }

  // Materialize the output in original input order.
  const out: DodgedLabel[] = new Array(peaks.length);
  for (let i = 0; i < indexed.length; i++) {
    const { idx, p, naturalPx } = indexed[i]!;
    const px = dodgedPx[i]!;
    const qLabel = px === naturalPx ? p.q : fromPx(px);
    out[idx] = {
      qPeak: p.q,
      yPeak: p.y,
      qLabel,
      yLabel: p.y - LABEL_OFFSET_PX,
      label: p.q.toPrecision(3),
      peakId: p.peakId,
    };
  }
  return out;
}
