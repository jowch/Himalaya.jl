/**
 * NeedsReviewBadge — review-mode header badge for stale comparisons
 * (Plan §Phase 9, Task 9.6; spec §Derived analysis state and staleness).
 *
 * Visual language: amber/warning tone (same family as `StaleIndicesBanner`),
 * NOT red/error. The badge is informational — the plot continues to render
 * the snapshotted state correctly; "stale" only means the underlying
 * analysis has drifted since submit time.
 *
 * Author-only clickable: clicking navigates the author to edit mode where
 * `loadDraftFromComparison` recomputes snapshots against the current cache.
 * Non-authors see the badge as informational; they can fork from the same
 * page if they want a refreshed version.
 *
 * Spec selectors:
 *   - container `data-testid="comparison-needs-review"`
 *   - `data-clickable={"true"|"false"}` reflects author gating
 */
import { useNavigate } from "react-router-dom";
import { useCurrentUserId } from "../hooks/useCurrentUserId";

export interface NeedsReviewBadgeProps {
  comparisonId: number;
  experimentId: number | undefined;
  /** Numeric user id of the author. Pass null for unauthored comparisons. */
  authorUserId: number | null;
}

const TOOLTIP_TEXT =
  "Underlying analysis has changed since this comparison was last submitted. Edit and re-submit to refresh.";

export function NeedsReviewBadge(
  { comparisonId, experimentId, authorUserId }: NeedsReviewBadgeProps,
): JSX.Element {
  const navigate = useNavigate();
  const currentUserId = useCurrentUserId();
  const isAuthor =
    authorUserId !== null
    && currentUserId !== undefined
    && currentUserId === authorUserId;

  const onClick = (): void => {
    if (!isAuthor) return;
    if (experimentId !== undefined) {
      navigate(`/experiments/${experimentId}/compare/${comparisonId}/edit`);
    } else {
      // Global scope still has an edit route via the all-listing page;
      // route shape mirrors `ComparePageEdit`'s navigation helpers.
      navigate(`/compare/all`);
    }
  };

  // Both clickable + non-clickable variants share visual language; only the
  // affordance + handler differ. The button form gives focus + Enter/Space
  // handling for free; the span form is explicitly aria-disabled so
  // screen readers don't announce a misleading interaction.
  if (isAuthor) {
    return (
      <button
        type="button"
        data-testid="comparison-needs-review"
        data-clickable="true"
        title={TOOLTIP_TEXT}
        onClick={onClick}
        className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full
                   border border-warning text-warning text-xs font-medium
                   hover:bg-bg-elevated cursor-pointer"
      >
        <span aria-hidden="true">⚠</span>
        <span>Needs Review</span>
      </button>
    );
  }
  return (
    <span
      data-testid="comparison-needs-review"
      data-clickable="false"
      role="status"
      aria-disabled="true"
      title={TOOLTIP_TEXT}
      className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full
                 border border-warning text-warning text-xs font-medium
                 cursor-default opacity-90"
    >
      <span aria-hidden="true">⚠</span>
      <span>Needs Review</span>
    </span>
  );
}
