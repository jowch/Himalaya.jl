import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "primary" | "ghost" | "danger";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
}

const variantClass: Record<ButtonVariant, string> = {
  primary:
    "bg-accent border border-accent text-white hover:brightness-110 " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghost:
    "text-fg-muted hover:text-fg hover:bg-bg-hover border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
  danger:
    "text-fg-muted hover:text-error border border-transparent " +
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
      className={`rounded-md px-2.5 py-1 transition-colors ${variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
}
