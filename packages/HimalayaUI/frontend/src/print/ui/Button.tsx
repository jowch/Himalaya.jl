import { forwardRef } from "react";
import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "success" | "ghost" | "danger" | "outline" | "outlineDanger" | "outlineSuccess" | "ghostInverse";

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
  // Constructive primary (the Keep verb): sage-filled sibling of `accent`,
  // legible on light AND dark (bg-ink) surfaces. Status semantics only —
  // sage is the design system's success hue (DESIGN.md T-4), never decoration.
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
  // Coloured-outline cull verbs (dock Drop/Keep, spec §5): status hue on the
  // border + text over a light plate, filling on hover.
  outlineDanger:
    "border border-error bg-plate text-error hover:bg-error hover:text-paper " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  outlineSuccess:
    "border border-success bg-plate text-success hover:bg-success hover:text-paper " +
    "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
  ghostInverse:
    "bg-transparent text-paper/70 hover:text-paper border border-transparent " +
    "focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
};

// Armed overrides the resting look with the terracotta active fill (mockup
// .tool-btn.armed). Kept as an override layer so any variant can be armed.
const armedClass =
  "bg-accent border border-accent text-paper hover:brightness-110";

// forwardRef so owners can hold the DOM button — e.g. LoupeSidePanel hands its
// "Manage" trigger to ManageTagsModal for focus-restore-on-close (mirrors how
// IconButton forwards its ref).
export const Button = forwardRef<HTMLButtonElement, ButtonProps>(function Button(
  {
    variant = "ghost",
    armed,
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
      // FO-RESCORE2 F12: a toggle button (one that passes `armed`) must keep
      // aria-pressed in BOTH states — render "false" when off, not drop the
      // attribute (APG toggle pattern). `armed` undefined → a non-toggle button,
      // which carries no aria-pressed at all. So bind the raw value: undefined
      // omits the attribute, false renders "false", true renders "true".
      aria-pressed={armed}
      className={`text-meta font-semibold whitespace-nowrap rounded-md px-2.5 py-1 transition-colors disabled:opacity-45 disabled:cursor-not-allowed ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
});
