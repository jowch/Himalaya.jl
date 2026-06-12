import { forwardRef } from "react";
import type { HTMLAttributes } from "react";
import { CheckCircle } from "./CheckCircle";

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
  /** Override the roving tabindex. When provided, replaces the default
   *  `disabled ? -1 : 0` (used by the Samples roving grid so the checkbox is
   *  the cell's single tab stop only when its cell is active). Absent →
   *  unchanged. */
  tabIndex?: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * Checkbox — interactive tri-state checkbox primitive (src/print/ui/).
 *
 * The VISUAL is the round `CheckCircle` badge (reused, rendered `decorative`),
 * wrapped in an interactive `<span role="checkbox">` that owns all semantics:
 * `aria-checked` (false / true / mixed for indeterminate), `data-checked` /
 * `data-indeterminate` (stable E2E selectors), `tabIndex`, Space/Enter
 * activation, a disabled guard, and the accent `focus-visible` ring. The inner
 * CheckCircle is `aria-hidden` so the role/label live only on the wrapper.
 *
 * Indeterminate note: CheckCircle has no indeterminate visual, and no real
 * consumer currently passes a truthy `indeterminate` (SampleTableRow only
 * plumbs the prop through). The prop is retained for API/test compatibility —
 * outer `aria-checked="mixed"` + `data-indeterminate` are honoured — but the
 * rendered badge falls back to the unchecked CheckCircle look.
 *
 * Appearance lives entirely in CheckCircle (token-driven); the wrapper carries
 * only placement + interaction classes — no inline colour literals here.
 */
export const Checkbox = forwardRef<HTMLSpanElement, CheckboxProps>(function Checkbox(
  {
    checked = false,
    indeterminate = false,
    disabled = false,
    onChange,
    tabIndex,
    className,
    ...rest
  },
  ref,
): JSX.Element {
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
      ref={ref}
      role="checkbox"
      aria-checked={ariaChecked}
      aria-disabled={disabled || undefined}
      tabIndex={tabIndex ?? (disabled ? -1 : 0)}
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
        "inline-flex items-center justify-center flex-shrink-0 rounded-full cursor-pointer select-none",
        "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
        disabled && "opacity-40 cursor-not-allowed",
        className,
      )}
      {...rest}
    >
      <CheckCircle checked={isChecked} decorative />
    </span>
  );
});
