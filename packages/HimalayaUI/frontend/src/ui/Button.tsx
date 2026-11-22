import { forwardRef } from "react";
import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "success" | "ghost" | "danger" | "outline" | "outlineAccent" | "outlineSuccess" | "ghostInverse";

interface ButtonProps extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "title"> {
  variant?: ButtonVariant;
  armed?: boolean;
  size?: ButtonSize;
}
export type ButtonSize = "md" | "lg";
const sizeClass: Record<ButtonSize, string> = {
  md: "text-meta px-2.5 py-1",
  lg: "text-body px-4 py-2",
};
const variantClass: Record<ButtonVariant, string> = {
  solid:
    "bg-ink border border-ink text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  accent:
    "bg-accent border border-accent text-paper hover:brightness-110 " +
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
    "border border-accent-soft bg-transparent text-accent hover:bg-accent hover:text-paper hover:border-accent " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outlineSuccess:
    "border border-success-soft bg-transparent text-success hover:bg-success hover:text-paper hover:border-success " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghostInverse:
    "bg-transparent text-paper/70 hover:text-paper border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
};
const armedClass =
  "bg-accent border border-accent text-paper hover:brightness-110";
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
      data-armed={armed ? "true" : undefined}
      aria-pressed={armed}
      className={`${sizeClass[size]} font-semibold whitespace-nowrap rounded-md transition-colors disabled:opacity-45 disabled:cursor-not-allowed ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
});