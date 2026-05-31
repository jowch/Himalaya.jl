function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface ScreenedMarkProps {
  done: boolean;
  className?: string;
}

export function ScreenedMark({ done, className }: ScreenedMarkProps): JSX.Element {
  return (
    <span
      data-testid="screened-mark"
      role="img"
      aria-label={done ? "Screened" : "Not screened"}
      data-done={done ? "true" : undefined}
      // 13px is a fixed semantic badge dimension (not body text), so the arbitrary
      // size is intentional and allowed in print/ui. The check stroke uses the paper
      // token (var(--color-paper)), never a raw hex.
      className={cx(
        "inline-block w-[13px] h-[13px] rounded-full flex-shrink-0 border",
        done ? "bg-ink border-ink" : "border-hair-strong",
        className,
      )}
    >
      {done && (
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
