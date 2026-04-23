import { useIndices, useReanalyzeExposure } from "../queries";

export interface StaleIndicesBannerProps {
  exposureId: number | undefined;
}

export function StaleIndicesBanner(
  { exposureId }: StaleIndicesBannerProps,
): JSX.Element | null {
  const q = useIndices(exposureId);
  // Always call the mutation hook so hook order is stable across renders.
  // When exposureId is undefined we never click the button anyway.
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
      <button
        className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
        disabled={reanalyze.isPending}
        onClick={() => reanalyze.mutate()}
      >
        {reanalyze.isPending ? "Re-analyzing…" : "Re-analyze"}
      </button>
    </div>
  );
}
