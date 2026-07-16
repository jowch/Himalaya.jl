/**
 * cursor — pure snap math for the WaterfallChart peak-q cursor. A pointer x
 * (container-relative px) snaps to the nearest detected peak's q across all
 * rows, within a px tolerance, reusing the plot engine's hitTestPeaks (which
 * ignores optimistic id<0 peaks and resolves ties to the later peak).
 */
import { hitTestPeaks } from "../plot/interaction";
import type { WaterfallRow } from "./waterfallModel";

/** All anchors across all rows as hit-test candidates. id+q suffice; equal q → equal px. */
export function cursorCandidates(rows: WaterfallRow[]): { id: number; q: number }[] {
  return rows.flatMap((r) => r.anchors.map((a) => ({ id: a.id, q: a.q })));
}

/** Snap a container-relative px to the nearest peak q within tolPx; null if none. */
export function snapToPeakQ(
  px: number,
  rows: WaterfallRow[],
  toPx: (q: number) => number,
  tolPx: number,
): number | null {
  const hit = hitTestPeaks(cursorCandidates(rows), px, toPx, tolPx);
  return hit ? hit.q : null;
}
