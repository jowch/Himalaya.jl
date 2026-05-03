import { useEffect, useState } from "react";
import { useIndices, useExposure, useReanalyzeExposure } from "../queries";
import { useExposureHasPendingPeakOps, useQueueOpStatus } from "../lib/queue/hooks";
import { Button } from "./ui";

// Stale state must persist this long before the banner appears. Suppresses
// flicker during cross-entity refetch races on a curation-confirmed event:
// the SSE post_state writes the new exposure hash + indices into the cache
// atomically, but observers fire on slightly different microtasks and we
// can briefly read mismatched (old, new) pairs. 150ms covers the React
// scheduling window without delaying genuine staleness reports noticeably.
const DEFAULT_STALE_DEBOUNCE_MS = 150;

export interface StaleIndicesBannerProps {
  exposureId: number | undefined;
  debounceMs?: number;
}

export function StaleIndicesBanner(
  { exposureId, debounceMs = DEFAULT_STALE_DEBOUNCE_MS }: StaleIndicesBannerProps,
): JSX.Element | null {
  const indicesQ = useIndices(exposureId);
  const exposureQ = useExposure(exposureId);
  const reanalyze = useReanalyzeExposure(exposureId ?? 0);
  // Drives the button's "Re-analyzing…" affordance. Reads from the
  // MutationCache by op kind + exposureId so the spinner survives across
  // hook re-renders and is consistent with other queue-driven UI.
  const reanalyzeStatus = useQueueOpStatus("reanalyze_exposure", exposureId);
  const reanalyzePending = reanalyzeStatus === "pending";
  // Hide the banner during any in-flight peak op for this exposure: the user
  // has just curated, the queue is mid-flight, and a brief hash mismatch is
  // expected. The op's onSuccess (or SSE post_state) will land the matching
  // hash and the banner stays hidden — no flicker.
  const hasPendingPeakOps = useExposureHasPendingPeakOps(exposureId);

  const indices = indicesQ.data ?? [];
  const exposure = exposureQ.data;
  const expectedHash = exposure?.analysis_inputs_hash ?? null;
  const stale = expectedHash
    ? indices.filter((i) => i.inputs_hash !== expectedHash)
    : [];
  const isStale = stale.length > 0;

  const [visible, setVisible] = useState(false);
  useEffect(() => {
    if (!isStale) {
      setVisible(false);
      return;
    }
    if (debounceMs <= 0) {
      setVisible(true);
      return;
    }
    const handle = setTimeout(() => setVisible(true), debounceMs);
    return () => clearTimeout(handle);
  }, [isStale, debounceMs, expectedHash]);

  if (exposureId === undefined) return null;
  if (!expectedHash) return null;
  if (hasPendingPeakOps) return null;
  if (!visible) return null;

  return (
    <div
      role="alert"
      className="flex items-center justify-between gap-4 px-3 py-2 mb-2 border border-warning text-fg bg-bg-elevated rounded-md"
    >
      <span>
        {stale.length} {stale.length === 1 ? "index is" : "indices are"} stale.
      </span>
      <Button
        variant="primary"
        disabled={reanalyzePending}
        onClick={() => reanalyze.mutate({})}
      >
        {reanalyzePending ? "Re-analyzing…" : "Re-analyze"}
      </Button>
    </div>
  );
}
