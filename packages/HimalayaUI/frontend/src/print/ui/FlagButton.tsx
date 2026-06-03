function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface FlagButtonProps {
  /** Parsed value, e.g. "1 : 0.25", rendered mono. */
  value: string;
  /** true → uncertain-parse look; false/omitted → resolved re-openable look. */
  flagged?: boolean;
  /** Toggles the read. */
  onClick?: () => void;
  /** PLACEMENT ONLY. */
  className?: string;
}

/**
 * The clickable parsed-value control on the series-scoping worksheet.
 *
 * - **flagged** — accent value with a dashed accent underline + a
 *   "▸ check the read" caption; clicking accepts the read.
 * - **ok** (default) — ink value with a dotted re-openable affordance that
 *   firms up on hover; clicking re-opens the value.
 */
export function FlagButton({ value, flagged, onClick, className }: FlagButtonProps): JSX.Element {
  return (
    <button
      type="button"
      data-testid="flag-button"
      data-state={flagged ? "flagged" : "ok"}
      onClick={onClick}
      {...(flagged ? {} : { title: "Click to re-open this value" })}
      className={cx(
        "group/fb block text-right font-mono cursor-pointer min-w-[92px] flex-shrink-0",
        flagged ? "text-accent" : "text-ink",
        className,
      )}
    >
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
        <span className="block text-[9px] font-bold uppercase tracking-wide text-accent mt-0.5">
          ▸ check the read
        </span>
      )}
    </button>
  );
}
