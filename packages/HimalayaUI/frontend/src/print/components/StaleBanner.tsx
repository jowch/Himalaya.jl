import { Button } from "../ui/Button";

export interface StaleBannerProps {
  /** Number of stale indices. <= 0 → renders nothing. */
  staleCount: number;
  /** When true the button is disabled and shows "Re-analyzing…". */
  pending?: boolean;
  /** Fired on button click (not called when pending). */
  onReanalyze: () => void;
  /** Placement-only className (margin, width, etc). */
  className?: string;
}

/**
 * Presentational banner shown when one or more indices are stale.
 *
 * Appearance is entirely token/primitive-based so lint:design passes.
 * The Focus page owns the debounce, cache reads, and pending-ops gating —
 * this component is pure presentation.
 *
 * NoticePill is not composed here: it has no warning tone and is a compact
 * pill — not a row-with-action container. A plain alert div with named token
 * classes is the correct choice.
 */
export function StaleBanner({
  staleCount,
  pending = false,
  onReanalyze,
  className = "",
}: StaleBannerProps): JSX.Element | null {
  if (staleCount <= 0) return null;

  const countLabel =
    staleCount === 1
      ? `${staleCount} index is stale.`
      : `${staleCount} indices are stale.`;

  return (
    <div
      role="alert"
      className={`flex items-center justify-between gap-4 px-3 py-2 border border-warning text-ink bg-plate rounded-md ${className}`.trim()}
    >
      <span>{countLabel}</span>
      <Button variant="solid" disabled={pending} onClick={onReanalyze}>
        {pending ? "Re-analyzing…" : "Re-analyze"}
      </Button>
    </div>
  );
}
