import type { IndexEntry } from "../api";
import { phaseColor } from "../phases";
import type { CombLine } from "../ui/WorkspacePlot";
import { storedFit } from "./indexing";

function fittedQs(c: IndexEntry): readonly number[] {
  const q1 = c.predicted_q[0] ?? 0;
  const fit = storedFit(c);
  if (fit.landed === 0 || !(q1 > 0) || !(fit.slope > 0)) return c.predicted_q;
  return c.predicted_q.map((q) => (q / q1) * fit.slope);
}

export function callCombs(
  candidates: readonly IndexEntry[],
  hoverIds: readonly number[] | undefined,
  assigned: ReadonlySet<number>,
  draft: { phase: string; qs: readonly number[] } | null,
): CombLine[] {
  if (draft) return [{ id: -1, qs: draft.qs, color: phaseColor(draft.phase) }];
  const hovered = candidates.filter((c) => hoverIds?.includes(c.id));
  if (hovered.length > 0) return hovered.map((c) => ({ id: c.id, qs: fittedQs(c), color: phaseColor(c.phase) }));
  return candidates.filter((c) => assigned.has(c.id)).map((c) => ({ id: c.id, qs: fittedQs(c), color: phaseColor(c.phase), faded: true }));
}
