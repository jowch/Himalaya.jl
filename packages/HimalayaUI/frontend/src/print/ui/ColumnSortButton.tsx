import { forwardRef } from "react";
import type { ButtonHTMLAttributes } from "react";

export type ColumnSortDir = "asc" | "desc";

export interface ColumnSortButtonProps
  extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "children"> {
  /** The visible column label, which is ALSO the button's accessible name. */
  label: string;
  /** "asc" | "desc" when this column is the active sort; undefined when it is
   *  sortable-but-inactive. Drives the caret glyph (the AT-facing state lives in
   *  the columnheader's `aria-sort`, set by the caller — NOT here). */
  active?: ColumnSortDir | undefined;
}

/**
 * The clickable label inside a sortable column header. Closed-look primitive:
 * appearance (soft-kicker label geometry/color, the resting/hover affordance,
 * the focus-visible ring, the direction caret) lives here; the consumer
 * (SheetTable) passes only the label, the active direction, and an onClick,
 * plus placement-only `className`.
 *
 * The caret (▲/▼) is `aria-hidden` decoration — it conveys nothing to AT, which
 * reads the column's state off the columnheader `aria-sort` the caller sets.
 * Inactive sortable headers show no caret at rest; the whole control gains a
 * subtle ink-darkening on hover/focus to read as actionable.
 *
 * The focus ring is `outline` (not a box-shadow ring) so it is never clipped by
 * the sticky header's own `overflow`/stacking — outline paints outside the box.
 *
 * Forwards its ref to the underlying `<button>` so the Samples roving grid can
 * register a sortable header as a roving cell and move focus to it.
 */
export const ColumnSortButton = forwardRef<HTMLButtonElement, ColumnSortButtonProps>(
  function ColumnSortButton({ label, active, className = "", ...props }, ref): JSX.Element {
  // ▲ ascending, ▼ descending. Only the active column shows a caret.
  const caret = active === "asc" ? "▲" : active === "desc" ? "▼" : null;

  return (
    <button
      ref={ref}
      type="button"
      data-sort-button="true"
      data-active={active ?? "none"}
      className={
        "group/sort flex w-full items-center gap-1 px-4 py-2.5 text-left " +
        "text-kicker text-kicker-soft transition-colors " +
        "hover:text-ink focus-visible:outline focus-visible:outline-2 " +
        "focus-visible:-outline-offset-2 focus-visible:outline-accent " +
        className
      }
      {...props}
    >
      <span>{label}</span>
      {/* caret slot: active → directional glyph; inactive → a faint neutral
          ▾ that only surfaces on hover/focus so the header reads as sortable. */}
      <span aria-hidden="true" className="text-ink-soft">
        {caret ?? (
          <span className="opacity-0 transition-opacity group-hover/sort:opacity-100 group-focus-visible/sort:opacity-100">
            ▾
          </span>
        )}
      </span>
    </button>
  );
});
