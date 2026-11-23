import type { IndexEntry } from "../api";
import { phaseColor } from "../phases";
import type { CombLine } from "../ui/WorkspacePlot";

export function callCombs(
  candidates: readonly IndexEntry[],
  hoverIds: readonly number[] | undefined,
  assigned: ReadonlySet<number>,
  draft: { phase: string; qs: readonly number[] } | null,
): CombLine[] {
  if (draft) return [{ id: -1, qs: draft.qs, color: phaseColor(draft.phase) }];
  const hovered = candidates.filter((c) => hoverIds?.includes(c.id));
  if (hovered.length > 0) return hovered.map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase) }));
  return candidates.filter((c) => assigned.has(c.id)).map((c) => ({ id: c.id, qs: c.predicted_q, color: phaseColor(c.phase), faded: true }));
}
