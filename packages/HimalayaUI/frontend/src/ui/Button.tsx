import { forwardRef } from "react";
import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "success" | "ghost" | "danger" | "outline" | "outlineAccent" | "outlineSuccess" | "ghostInverse";

interface ButtonProps extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "title"> {
  variant?: ButtonVariant;
  armed?: boolean;
  size?: ButtonSize;
}
export type ButtonSize = "sm" | "md" | "lg";
export const buttonSizeClass: Record<ButtonSize, string> = {
  sm: "inline-flex h-6 items-center justify-center gap-1 rounded-sm px-2 text-caption",
  md: "inline-flex h-7 items-center justify-center gap-1.5 rounded-md px-2.5 text-meta font-semibold",
  lg: "inline-flex h-9.5 items-center justify-center gap-2 rounded-md px-4 text-body font-semibold",
};
const variantClass: Record<ButtonVariant, string> = {
  solid:
    "bg-ink border border-ink text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  accent:
    "bg-accent-fill border border-accent-fill text-on-accent hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  success:
    "bg-success border border-success text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghost:
    "text-ink-soft hover:text-ink hover:bg-paper-sunk border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
  danger:
    "text-error border border-transparent hover:bg-error hover:text-paper hover:border-error " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outline:
    "border border-hair-strong bg-plate text-ink hover:bg-paper-sunk " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outlineAccent:
    "border border-accent-soft bg-transparent text-accent hover:bg-accent-fill hover:text-on-accent hover:border-accent-fill " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outlineSuccess:
    "border border-success-soft bg-transparent text-success hover:bg-success hover:text-paper hover:border-success " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghostInverse:
    "bg-transparent text-paper/70 hover:text-paper border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
};
const armedClass =
  "bg-accent-fill border border-accent-fill text-on-accent hover:brightness-110";
export const Button = forwardRef<HTMLButtonElement, ButtonProps>(function Button(
  {
    variant = "ghost",
    armed,
    size = "md",
    className = "",
    children,
    ...props
  },
  ref,
): JSX.Element {
  return (
    <button
      ref={ref}
      data-variant={variant}
      data-size={size}
      data-armed={armed ? "true" : undefined}
      aria-pressed={armed}
      className={`${buttonSizeClass[size]} whitespace-nowrap transition-colors disabled:opacity-45 disabled:cursor-not-allowed ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
});