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
 * - **ok** (default) — the read will be committed; ink value with a dotted
 *   affordance that firms up on hover. Clicking skips the read.
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
      aria-pressed={flagged === true}
      onClick={onClick}
      title={flagged ? "Restore this read" : "Skip this read"}
      className={cx(
        "group/fb block text-right font-mono cursor-pointer min-w-[92px] flex-shrink-0",
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
            : "border-dotted border-transparent group-hover/fb:border-hair-strong",
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
