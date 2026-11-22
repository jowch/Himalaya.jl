import { forwardRef } from "react";
import type { HTMLAttributes } from "react";
import { CheckCircle } from "./CheckCircle";
import { cx } from "../lib/cx";

export interface CheckboxProps extends Omit<HTMLAttributes<HTMLSpanElement>, "onChange" | "title"> {
  checked?: boolean;
  indeterminate?: boolean;
  disabled?: boolean;
  onChange?: (checked: boolean) => void;
  "aria-label"?: string;
  tabIndex?: number;
  className?: string;
}
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