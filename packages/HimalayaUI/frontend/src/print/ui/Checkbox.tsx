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
  /** Override tabIndex. When provided, replaces the default `disabled ? -1 : 0`.
   *  Absent → unchanged. */
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
      data-hit-area="24"
      onClick={activate}
      onKeyDown={(e) => {
        if (e.key === " " || e.key === "Enter") {
          e.preventDefault();
          activate();
        }
      }}
      className={cx(
        // p-1.5 grows the interactive HIT AREA to 13px disc + 12px padding = 25px
        // (>= the 24x24 WCAG 2.5.8 minimum) while -m-1.5 cancels it in layout so
        // the visible disc and the surrounding column geometry are unchanged. The
        // disc's own size lives in CheckCircle and is untouched.
        "inline-flex items-center justify-center flex-shrink-0 rounded-full cursor-pointer select-none",
        "p-1.5 -m-1.5",
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
