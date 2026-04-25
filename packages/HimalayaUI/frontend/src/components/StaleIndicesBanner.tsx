import { useIndices, useReanalyzeExposure } from "../queries";
import { Button } from "./ui";

export interface StaleIndicesBannerProps {
  exposureId: number | undefined;
}

export function StaleIndicesBanner(
  { exposureId }: StaleIndicesBannerProps,
): JSX.Element | null {
  const q = useIndices(exposureId);
  const reanalyze = useReanalyzeExposure(exposureId ?? 0);

  if (exposureId === undefined) return null;
  const indices = q.data ?? [];
  const stale = indices.filter((i) => i.status === "stale");
  if (stale.length === 0) return null;

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
