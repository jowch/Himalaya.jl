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
 * Dismissal: trigger toggle, Esc keydown, or outside-click (issue #75).
 */
import { useEffect, useRef, useState } from "react";
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
  const triggerRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);
  const forksQ = useComparisonForks(comparisonId);
  const forks: ComparisonSummary[] = forksQ.data ?? [];

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") setOpen(false);
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [open]);

  useEffect(() => {
    if (!open) return;
    const onDown = (e: MouseEvent): void => {
      const target = e.target as Node | null;
      if (target === null) return;
      if (panelRef.current?.contains(target)) return;
      if (triggerRef.current?.contains(target)) return;
      setOpen(false);
    };
    document.addEventListener("mousedown", onDown);
    return () => document.removeEventListener("mousedown", onDown);
  }, [open]);

  return (
    <div className="relative inline-block">
      <button
        ref={triggerRef}
        type="button"
        data-testid="comparison-forks-trigger"
        aria-expanded={open}
        aria-haspopup="true"
        onClick={() => setOpen((v) => !v)}
        className="px-2 py-0.5 rounded border border-border text-fg-muted text-xs
                   hover:bg-bg-elevated"
      >
        Forks ({forks.length}) →
      </button>
      {open && (
        <div
          ref={panelRef}
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
