function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface FlagButtonProps {
  /** Parsed value, e.g. "1 : 0.25", rendered mono. */
  value: string;
  /** true → the user skipped this read from the batch write; false/omitted → the read will be committed. */
  flagged?: boolean;
  /** Toggles whether the read is skipped. */
  onClick?: () => void;
  /** PLACEMENT ONLY. */
  className?: string;
}

/**
 * The clickable parsed-value control on the series-scoping worksheet. It is a
 * skip toggle: clicking excludes (or restores) this read from the batch write.
 *
 * - **flagged** — the user skipped this read from the batch write; the accent
 *   look (accent value, dashed underline, "▸ skipped" caption) marks the
 *   exclusion. Clicking restores the read.
 * - **ok** (default) — the read will be committed; ink value with a standing
 *   faint dotted underline that firms up on hover. Clicking skips the read.
 *
 * Discovery (SC-SKIPDISC): the underline is PERSISTENT, not hover-only — the
 * skip toggle must be findable without hovering blindly. The contract is
 * exposed as `data-affordance="persistent"` (jsdom can't compute Tailwind
 * styles, so tests pin the attribute, not the class). Keyboard users get the
 * house focus-visible accent outline (same tokens as the reorder grip).
 *
 * The `flagged` prop name and the `data-state="ok"/"flagged"` values are KEPT
 * as stable API/e2e keys (ScopeSampleRow's row-wash keys off the same name) —
 * only the user-facing vocabulary is skip/skipped. State is exposed via
 * `aria-pressed`; the accessible name ("Skip this read: {value}") stays
 * constant across states.
 */
export function FlagButton({ value, flagged, onClick, className }: FlagButtonProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid="flag-button"
      data-state={flagged ? "flagged" : "ok"}
      data-affordance="persistent"
      aria-pressed={flagged === true}
      onClick={onClick}
      title={flagged ? "Restore this read" : "Skip this read"}
      className={cx(
        "group/fb block text-right font-mono cursor-pointer min-w-[92px] flex-shrink-0",
        "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
        flagged ? "text-accent" : "text-ink",
        className,
      )}
    >
      <span className="sr-only">Skip this read: </span>
      <span
        className={cx(
          "text-data font-bold border-b pb-px",
          flagged
            ? "border-dashed border-accent"
            : "border-dotted border-hair group-hover/fb:border-hair-strong",
        )}
      >
        {value}
      </span>
      {flagged && (
        <span
          aria-hidden="true"
          className="block text-[9px] font-bold uppercase tracking-wide text-accent mt-0.5"
        >
          ▸ skipped
        </span>
      )}
    </button>
  );
}
