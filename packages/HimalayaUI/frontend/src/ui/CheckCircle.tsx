import { cx } from "../lib/cx";

interface CheckCircleProps {
  checked: boolean;
  label?: string;
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