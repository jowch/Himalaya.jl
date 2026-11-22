import { forwardRef } from "react";
import type { ButtonHTMLAttributes, ReactNode } from "react";
import { Tooltip } from "./Tooltip";
import { cx } from "../lib/cx";

export type IconButtonTone = "ghost" | "accent" | "danger";

export interface IconButtonProps
  extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "aria-label"> {
  label: string;
  tone?: IconButtonTone;
  boxed?: boolean;
  dismiss?: boolean;
  children?: ReactNode;
}
const toneClass: Record<IconButtonTone, string> = {
  ghost: "text-ink-faint hover:text-ink",
  accent: "text-ink-faint hover:text-print-accent",
  danger: "text-error hover:bg-error/10",
};
export const IconButton = forwardRef<HTMLButtonElement, IconButtonProps>(function IconButton(
  {
    label,
    tone = "ghost",
    boxed = false,
    dismiss = false,
    className = "",
    children,
    type = "button",
    title,
    ...props
  },
  ref,
): JSX.Element {
  return (
    <Tooltip content={title ?? label}>
      <button
        ref={ref}
        type={type}
        aria-label={label}
        data-tone={tone}
        className={cx(
          "icon-button relative inline-flex items-center justify-center",
          "leading-none transition-colors",
          boxed
            ? "h-[26px] w-[26px] rounded-md border border-hair-strong bg-paper"
            : "rounded p-1",
          "focus-visible:outline focus-visible:outline-2",
          "focus-visible:outline-offset-2 focus-visible:outline-accent",
          "disabled:cursor-not-allowed disabled:opacity-30",
          toneClass[tone],
          className,
        )}
        {...props}
      >
        {dismiss ? "×" : children}
      </button>
    </Tooltip>
  );
});