function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface CheckCircleProps {
  checked: boolean;
  /** Accessible name override. Defaults to "Selected"/"Not selected"; a
   *  screened-status use can pass label={done ? "Screened" : "Not screened"}. */
  label?: string;
  /** When true, the badge is purely decorative: no `role`/`aria-label`, and
   *  `aria-hidden`. Used when an interactive wrapper (e.g. the Checkbox
   *  primitive) already provides the checkbox semantics — avoids doubling up
   *  roles/labels on nested elements. */
  decorative?: boolean;
  className?: string;
}

export function CheckCircle({ checked, label, decorative, className }: CheckCircleProps): JSX.Element {
  return (
    <span
      data-testid="check-circle"
      role={decorative ? undefined : "img"}
      aria-label={decorative ? undefined : (label ?? (checked ? "Selected" : "Not selected"))}
      aria-hidden={decorative ? "true" : undefined}
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
