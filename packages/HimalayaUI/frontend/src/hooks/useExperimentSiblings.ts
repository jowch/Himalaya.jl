import { useMemo } from "react";
import { useCorpusSamples } from "../queries";
import { useAppState } from "../state";
import {
  deriveExperimentSiblings,
  type ExperimentSiblings,
} from "../lib/sample/experimentSiblings";

/**
 * useExperimentSiblings — the live wiring of `deriveExperimentSiblings`:
 * corpus cache (TanStack) × active sample (Zustand). F5: shared by the topbar
 * sample stepper AND the Focus sample step (FocusPage, via the page's interaction
 * declaration) so both step through the identical sibling list. Deliberately does NOT depend
 * on `activeExperimentId` (only the NavModal picker / recoverFromStaleUrl write
 * that) — the experiment scope comes from the active sample's own `experiment_id`.
 *
 * Memoized so effect consumers re-subscribe only when the corpus list or the
 * active sample actually changes.
 */
export function useExperimentSiblings(): ExperimentSiblings {
  const corpusQ = useCorpusSamples();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  return useMemo(
    () => deriveExperimentSiblings(corpusQ.data, activeSampleId),
    [corpusQ.data, activeSampleId],
  );
}
