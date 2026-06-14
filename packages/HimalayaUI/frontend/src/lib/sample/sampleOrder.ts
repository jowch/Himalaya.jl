import type { CorpusSample } from "../../api";

/**
 * The ordered sample-id list the loupe steps through, shared by LoupePage (its
 * keyboard `[`/`]` + gotoSample) and CorpusTopbar (the sample stepper), so the
 * two can never disagree about prev/next.
 *
 * Priority: the order carried in router `state.sampleOrder` (the contact-sheet
 * order the loupe was opened with, preserved across steps) when it still
 * contains the active sample; otherwise the corpus scoped to the active beamtime
 * (experiment), in corpus order.
 */
export function resolveSampleOrder(
  corpus: CorpusSample[],
  beamtime: number | undefined,
  sampleId: number,
  sampleOrderState: number[] | undefined,
): number[] {
  if (sampleOrderState && sampleOrderState.includes(sampleId)) {
    return sampleOrderState;
  }
  const scoped =
    beamtime === undefined ? corpus : corpus.filter((s) => s.experiment_id === beamtime);
  return scoped.map((s) => s.id);
}

/** Prev/next ids + the 0-based index for `sampleId` within `ordered`. */
export function sampleNeighbors(
  ordered: number[],
  sampleId: number,
): { index: number; prevId: number | undefined; nextId: number | undefined } {
  const index = ordered.indexOf(sampleId);
  return {
    index,
    prevId: index > 0 ? ordered[index - 1] : undefined,
    nextId: index >= 0 && index < ordered.length - 1 ? ordered[index + 1] : undefined,
  };
}
