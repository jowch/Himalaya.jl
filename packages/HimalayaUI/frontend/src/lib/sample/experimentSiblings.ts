// experimentSiblings.ts — single source of truth for "which samples does the
// focus surface step between". F5: the app shell stepper and the
// Focus `[`/`]` sample step (FocusPage, via useShortcuts) both consume this
// derivation, so the two can never disagree (alignment invariant — the same
// rule as the shared SAMPLE_TABLE_COLS grid constant).
import type { CorpusSample } from "../../api";

export interface ExperimentSiblings {
  /** The corpus row for the active sample; undefined when the cache is cold
   *  or the id is unknown — every other field degrades to empty with it. */
  activeSample: CorpusSample | undefined;
  /** The active sample's experiment-siblings, in corpus order. */
  siblings: CorpusSample[];
  /** Index of the active sample within `siblings`; -1 when unresolved. */
  index: number;
  /** Previous/next sibling. No wrap: undefined at the ends — and undefined
   *  whenever the derivation is unresolved, so consumers no-op gracefully. */
  prev: CorpusSample | undefined;
  next: CorpusSample | undefined;
}

/**
 * Derive the active sample's experiment-siblings from the corpus list:
 * resolve the active sample, filter the corpus to its experiment (preserving
 * corpus order), and locate prev/next without wrapping.
 *
 * Cold cache (`samples === undefined`) or an unknown/absent `activeSampleId`
 * yields the empty derivation — callers need no extra guards.
 */
export function deriveExperimentSiblings(
  samples: CorpusSample[] | undefined,
  activeSampleId: number | undefined,
): ExperimentSiblings {
  const activeSample =
    activeSampleId !== undefined
      ? samples?.find((s) => s.id === activeSampleId)
      : undefined;
  const siblings =
    activeSample !== undefined
      ? (samples ?? []).filter((s) => s.experiment_id === activeSample.experiment_id)
      : [];
  const index =
    activeSample !== undefined
      ? siblings.findIndex((s) => s.id === activeSample.id)
      : -1;
  const prev = index > 0 ? siblings[index - 1] : undefined;
  const next =
    index >= 0 && index < siblings.length - 1 ? siblings[index + 1] : undefined;
  return { activeSample, siblings, index, prev, next };
}
