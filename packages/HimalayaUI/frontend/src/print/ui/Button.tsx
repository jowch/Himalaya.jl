import { forwardRef } from "react";
import type { ButtonHTMLAttributes } from "react";

export type ButtonVariant = "solid" | "accent" | "success" | "ghost" | "danger" | "outline" | "outlineAccent" | "outlineSuccess" | "ghostInverse";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: ButtonVariant;
  /** Toggled/active tool-button state (focus-plot "+ Peak" armed): terracotta
   *  fill, paper text, `aria-pressed`. Distinct from `variant="accent"` (a
   *  primary action, not a toggle). */
  armed?: boolean;
  /** "md" (default): the dense chrome button. "lg": a prominent primary CTA
   *  (~50% larger box + body-size label) for empty-state / first-run calls to
   *  action, where the button is the focal point rather than dense chrome. */
  size?: ButtonSize;
}

export type ButtonSize = "md" | "lg";

// Size owns the box geometry + label scale; variant owns colour. Kept orthogonal
// so any variant can be sized up for an empty-state hero CTA.
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
  // Coloured-outline cull verbs (dock Drop/Keep, spec §5): the verb's hue
  // (accent for Drop, success for Keep) on the border + text over a light
  // plate, filling on hover. Accent Drop matches CullBar (one Drop hue app-wide).
  outlineAccent:
    "border border-accent bg-plate text-accent hover:bg-accent hover:text-paper " +
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
      // FO-RESCORE2 F12: a toggle button (one that passes `armed`) must keep
      // aria-pressed in BOTH states — render "false" when off, not drop the
      // attribute (APG toggle pattern). `armed` undefined → a non-toggle button,
      // which carries no aria-pressed at all. So bind the raw value: undefined
      // omits the attribute, false renders "false", true renders "true".
      aria-pressed={armed}
      className={`${sizeClass[size]} font-semibold whitespace-nowrap rounded-md transition-colors disabled:opacity-45 disabled:cursor-not-allowed ${armed ? armedClass : variantClass[variant]} ${className}`}
      {...props}
    >
      {children}
    </button>
  );
});
