import type { ReactNode } from "react";

interface KbKeyProps {
  children: ReactNode;
  /** `frost` — translucent currentColor fill for a chip sitting ON a filled
   *  colored button (the §7 Focus-primary `↵`): it inherits the button's text
   *  color rather than the plate look. Default `plate`. */
  variant?: "plate" | "frost";
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

// Per-variant appearance. `plate` is the standalone chip (its own surface);
// `frost` is a translucent currentColor fill so the chip reads as part of a
// filled colored button and inherits that button's text color (§7 frost).
const VARIANT: Record<NonNullable<KbKeyProps["variant"]>, string> = {
  plate: "bg-plate text-ink-soft border-hair-strong",
  frost: "bg-current/15 text-current border-current/40",
};

/** Keyboard key cap — a semantic <kbd> badge for shortcut hints (e.g. ⌘K).
 *  Mono is the right voice for a key cap (a measured/literal token, not prose).
 *  Static label: no interactive state, so checklist C (interaction states) and
 *  D (touch target) are N/A here.
 *
 *  The `border-b-2` is a UNIFORM border-colored border that is simply thicker
 *  (2px) at the bottom — the physical-key 3D affordance of a raised key cap.
 *  It is NOT a colored accent side-stripe (checklist G ban): the color is
 *  identical to the rest of the border, only the bottom edge weight differs. */
export function KbKey({ children, variant = "plate", className = "" }: KbKeyProps): JSX.Element {
  return (
    <kbd
      data-testid="kbkey"
      className={cx(
        "inline-block font-mono text-xs rounded-sm px-1.5 py-px",
        VARIANT[variant],
        // Uniform border, just thicker (2px) at the bottom edge — the raised-key
        // 3D affordance, NOT a colored accent stripe (checklist G).
        "border border-b-2",
        className,
      )}
    >
      {children}
    </kbd>
  );
}
