function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface CheckCircleProps {
  checked: boolean;
  /** Accessible name override. Defaults to "Selected"/"Not selected"; a
   *  screened-status use can pass label={done ? "Screened" : "Not screened"}. */
  label?: string;
  className?: string;
}

export function CheckCircle({ checked, label, className }: CheckCircleProps): JSX.Element {
  return (
    <span
      data-testid="check-circle"
      role="img"
      aria-label={label ?? (checked ? "Selected" : "Not selected")}
      data-checked={checked ? "true" : undefined}
      // 13px is a fixed semantic badge dimension (not body text), so the arbitrary
      // size is intentional and allowed in print/ui. The check stroke uses the paper
      // token (var(--color-paper)), never a raw hex.
      className={cx(
        "inline-block w-[13px] h-[13px] rounded-full flex-shrink-0 border",
        checked ? "bg-ink border-ink" : "border-hair-strong",
        className,
      )}
    >
      {checked && (
        <svg viewBox="0 0 13 13" className="w-full h-full" aria-hidden="true">
          <path
            d="M3.4 6.8l2.1 2.1 4.2-4.6"
            fill="none"
            stroke="var(--color-paper)"
            strokeWidth={1.7}
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      )}
    </span>
  );
}
