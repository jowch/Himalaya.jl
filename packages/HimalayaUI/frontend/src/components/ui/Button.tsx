import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "ghost" | "danger";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
}

const variantClass: Record<ButtonVariant, string> = {
  solid:
    "bg-ink border border-ink text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  accent:
    "bg-accent border border-accent text-paper hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghost:
    "text-ink-soft hover:text-ink hover:bg-paper-sunk border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
  danger:
    "text-ink-soft hover:text-error border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
};

export function Button({
  variant = "ghost",
  className = "",
  children,
  ...props
}: ButtonProps): JSX.Element {
  return (
    <button
      data-variant={variant}
      className={`rounded-md px-2.5 py-1 transition-colors ${variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
}
