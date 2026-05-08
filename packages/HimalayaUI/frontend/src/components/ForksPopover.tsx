/**
 * ForksPopover — "Forks (N) →" trigger + dropdown listing the comparison's
 * child forks (Plan §Phase 11, Task 11.3; spec §Authorship and forking).
 *
 * Click-to-toggle popover. Lists each fork with its title + author + a link
 * to that fork's review page. Empty state surfaces "No forks yet".
 *
 * Selectors:
 *   - `data-testid="comparison-forks-trigger"` on the button.
 *   - `data-testid="comparison-forks-popover"` on the dropdown panel.
 *   - `data-testid="comparison-forks-row"` on each fork row.
 *
 * The popover is intentionally lightweight (no focus trap; not a modal).
 * Clicking the trigger again closes it. A future enhancement could add
 * outside-click dismiss, but the trigger toggle is enough for v1.
 */
import { useState } from "react";
import { Link } from "react-router-dom";
import { useComparisonForks } from "../queries";
import { comparePath, type CompareScope } from "../lib/comparison/routes";
import type { ComparisonSummary } from "../api";

export interface ForksPopoverProps {
  comparisonId: number;
  /** Experiment context (undefined → global /compare/all routing for forks). */
  experimentId: number | undefined;
  /**
   * Explicit scope override; when omitted, falls back to eid-derived. The
   * Compare page passes scope so /compare/all/:id's fork list deep-links
   * each child to /compare/all/:childId rather than possibly jumping to
   * an experiment route.
   */
  scope?: CompareScope;
}

export function ForksPopover({
  comparisonId, experimentId, scope,
}: ForksPopoverProps): JSX.Element {
  const resolvedScope: CompareScope = scope
    ?? (experimentId !== undefined ? "experiment" : "all");
  const [open, setOpen] = useState(false);
  const forksQ = useComparisonForks(comparisonId);
  const forks: ComparisonSummary[] = forksQ.data ?? [];

  return (
    <div className="relative inline-block">
      <button
        type="button"
        data-testid="comparison-forks-trigger"
        onClick={() => setOpen((v) => !v)}
        className="px-1.5 py-0.5 rounded text-xs text-fg-dim hover:text-fg
                   hover:bg-bg-hover border border-transparent hover:border-border"
      >
        Forks ({forks.length}) →
      </button>
      {open && (
        <div
          data-testid="comparison-forks-popover"
          className="absolute z-50 mt-1 right-0 min-w-[260px] max-h-[280px]
                     overflow-y-auto card border border-border bg-bg-elevated
                     shadow-lg p-2"
        >
          {forksQ.isLoading ? (
            <div
              data-testid="comparison-forks-loading"
              className="text-fg-muted text-xs px-2 py-1 italic"
            >
              Loading forks…
            </div>
          ) : forks.length === 0 ? (
            <div className="text-fg-muted text-xs px-2 py-1">No forks yet</div>
          ) : (
            <ul className="flex flex-col gap-1">
              {forks.map((f) => (
                <li
                  key={f.id}
                  data-testid="comparison-forks-row"
                  className="text-xs"
                >
                  <Link
                    to={comparePath({
                      scope: resolvedScope,
                      eid: experimentId,
                      id: f.id,
                    })}
                    className="block px-2 py-1 rounded hover:bg-bg
                               text-fg hover:text-accent"
                    onClick={() => setOpen(false)}
                  >
                    <div className="font-medium truncate">
                      {f.title || `Comparison #${f.id}`}
                    </div>
                  </Link>
                </li>
              ))}
            </ul>
          )}
        </div>
      )}
    </div>
  );
}
