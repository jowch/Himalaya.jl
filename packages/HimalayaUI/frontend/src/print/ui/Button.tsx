import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "ghost" | "danger" | "outline" | "ghostInverse";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
  /** Toggled/active tool-button state (focus-plot "+ Peak" armed): terracotta
   *  fill, paper text, `aria-pressed`. Distinct from `variant="accent"` (a
   *  primary action, not a toggle). */
  armed?: boolean;
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
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
  danger:
    "text-error border border-transparent hover:bg-error hover:text-paper hover:border-error " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outline:
    "border border-hair-strong bg-plate text-ink hover:bg-paper-sunk " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghostInverse:
    "bg-transparent text-paper/70 hover:text-paper border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
};

// Armed overrides the resting look with the terracotta active fill (mockup
// .tool-btn.armed). Kept as an override layer so any variant can be armed.
const armedClass =
  "bg-accent border border-accent text-paper hover:brightness-110";

export function Button({
  variant = "ghost",
  armed = false,
  className = "",
  children,
  ...props
}: ButtonProps): JSX.Element {
  return (
    <button
      data-variant={variant}
      data-armed={armed ? "true" : undefined}
      aria-pressed={armed ? true : undefined}
      className={`text-meta font-semibold rounded-md px-2.5 py-1 transition-colors disabled:opacity-45 disabled:cursor-not-allowed ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
}
