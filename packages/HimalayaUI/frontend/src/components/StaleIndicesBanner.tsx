import { useEffect, useState } from "react";
import { useIndices, useExposure, useReanalyzeExposure } from "../queries";
import { Button } from "./ui";

// Stale state must persist this long before the banner appears. Suppresses
// the flicker during the normal `peak op → autoReanalyze → invalidate` chain
// while still surfacing genuine failures (silent autoReanalyze errors,
// remote-actor edits where the actor's tab closes mid-edit).
const DEFAULT_STALE_DEBOUNCE_MS = 2000;

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
        disabled={reanalyze.isPending}
        onClick={() => reanalyze.mutate()}
      >
        {reanalyze.isPending ? "Re-analyzing…" : "Re-analyze"}
      </Button>
    </div>
  );
}
