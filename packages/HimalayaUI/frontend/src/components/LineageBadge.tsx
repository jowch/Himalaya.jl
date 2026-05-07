/**
 * LineageBadge — surfaces fork lineage in the review-mode header
 * (Plan §Phase 11, Task 11.3; spec §Authorship and forking).
 *
 * Renders nothing for non-fork comparisons. When the comparison is a fork
 * we surface either:
 *   - "Forked from <parent title> [view parent →]" when the parent still
 *     exists (`forked_from_id` is set).
 *   - "Forked from a deleted comparison" when the parent has been deleted
 *     (`forked_from_id IS NULL` post-cascade, but `forked_at_hash` is set
 *     because that field has no FK to clear).
 *
 * `forked_from_title` is denormalized into the GET payload by
 * `fetch_comparison_with_members` for cheap rendering without a second
 * fetch — falls back to "comparison" if missing.
 *
 * Selectors mirror spec §Selectors:
 *   - `data-testid="comparison-lineage"`
 *   - `data-parent-id={id|"deleted"}`
 *   - `data-testid="comparison-lineage-view-parent"` on the link (only
 *     when the parent still exists).
 */
import { Link } from "react-router-dom";
import type { Comparison } from "../api";

export interface LineageBadgeProps {
  comparison: Comparison;
  /** Experiment context (undefined → global /compare/all). */
  experimentId: number | undefined;
}

export function LineageBadge({
  comparison, experimentId,
}: LineageBadgeProps): JSX.Element | null {
  const { forked_from_id, forked_at_hash, forked_from_title } = comparison;

  // Not a fork at all.
  if (forked_from_id === null && forked_at_hash === null) return null;

  // Parent has been deleted: the FK ON DELETE SET NULL clears the id but
  // we still have the hash. Surface a degraded badge with no link.
  if (forked_from_id === null) {
    return (
      <span
        data-testid="comparison-lineage"
        data-parent-id="deleted"
        className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full
                   border border-border text-fg-muted text-xs"
        title="The parent comparison has been deleted."
      >
        Forked from a deleted comparison
      </span>
    );
  }

  // Live parent: build the link and surface "Forked from <title>".
  const parentHref = experimentId !== undefined
    ? `/experiments/${experimentId}/compare/${forked_from_id}`
    : "/compare/all";
  const parentTitle = forked_from_title ?? "comparison";
  return (
    <span
      data-testid="comparison-lineage"
      data-parent-id={String(forked_from_id)}
      className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full
                 border border-border text-fg-muted text-xs"
    >
      <span>
        Forked from <span className="text-fg">{parentTitle}</span>
      </span>
      <Link
        data-testid="comparison-lineage-view-parent"
        to={parentHref}
        className="text-accent hover:underline ml-1"
      >
        view parent →
      </Link>
    </span>
  );
}
