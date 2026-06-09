import type { HTMLAttributes } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface CheckboxProps extends Omit<HTMLAttributes<HTMLSpanElement>, "onChange"> {
  checked?: boolean;
  indeterminate?: boolean;
  disabled?: boolean;
  onChange?: (checked: boolean) => void;
  /** Accessible name — required when no visible label wraps this. */
  "aria-label"?: string;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * Checkbox — a closed-look tri-state checkbox primitive (src/print/ui/).
 *
 * Implemented as a `<span role="checkbox">` rather than a native `<input
 * type="checkbox">` so the visual ring is fully token-driven (bg-ink,
 * border-hair-strong, accent) and no browser reset is needed. The span
 * receives `tabIndex={0}`, responds to Space/Enter, and announces via
 * `aria-checked` (false / true / mixed for indeterminate). `data-checked`
 * and `data-indeterminate` are stable E2E selectors.
 *
 * Appearance is intentionally minimal: a 16px box, border-hair-strong resting
 * ring, bg-ink fill + paper SVG check when checked, accent fill + dash when
 * indeterminate. No inline colour literals or arbitrary sizes — all from @theme.
 */
export function Checkbox({
  checked = false,
  indeterminate = false,
  disabled = false,
  onChange,
  className,
  ...rest
}: CheckboxProps): JSX.Element {
  const isChecked = !indeterminate && checked;
  const ariaChecked: "true" | "false" | "mixed" = indeterminate
    ? "mixed"
    : checked
      ? "true"
      : "false";

  function activate(): void {
    if (!disabled) onChange?.(!checked);
  }

  return (
    <span
      role="checkbox"
      aria-checked={ariaChecked}
      aria-disabled={disabled || undefined}
      tabIndex={disabled ? -1 : 0}
      data-checked={isChecked ? "true" : "false"}
      data-indeterminate={indeterminate ? "true" : undefined}
      onClick={activate}
      onKeyDown={(e) => {
        if (e.key === " " || e.key === "Enter") {
          e.preventDefault();
          activate();
        }
      }}
      className={cx(
        "inline-flex items-center justify-center flex-shrink-0 w-4 h-4 rounded border cursor-pointer select-none transition-colors",
        "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
        disabled && "opacity-40 cursor-not-allowed",
        isChecked
          ? "bg-ink border-ink"
          : indeterminate
            ? "bg-accent border-accent"
            : "bg-plate border-hair-strong hover:border-ink",
        className,
      )}
      {...rest}
    >
      {isChecked && (
        <svg viewBox="0 0 16 16" className="w-full h-full" aria-hidden="true">
          <path
            d="M4 8.2l2.6 2.6 5.4-5.6"
            fill="none"
            stroke="var(--color-paper)"
            strokeWidth={1.8}
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      )}
      {indeterminate && (
        <svg viewBox="0 0 16 16" className="w-full h-full" aria-hidden="true">
          <line
            x1="4"
            y1="8"
            x2="12"
            y2="8"
            stroke="var(--color-paper)"
            strokeWidth={2}
            strokeLinecap="round"
          />
        </svg>
      )}
    </span>
  );
}
