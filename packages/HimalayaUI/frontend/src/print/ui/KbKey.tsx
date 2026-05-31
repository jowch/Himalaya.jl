import type { ReactNode } from "react";

interface KbKeyProps {
  children: ReactNode;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Keyboard key cap — a semantic <kbd> badge for shortcut hints (e.g. ⌘K).
 *  Mono is the right voice for a key cap (a measured/literal token, not prose).
 *  Static label: no interactive state, so checklist C (interaction states) and
 *  D (touch target) are N/A here.
 *
 *  The `border-b-2` is a UNIFORM border-hair-strong-colored border that is simply
 *  thicker (2px) at the bottom — the physical-key 3D affordance of a raised key cap.
 *  It is NOT a colored accent side-stripe (checklist G ban): the color is identical
 *  to the rest of the border, only the bottom edge weight differs. */
export function KbKey({ children, className = "" }: KbKeyProps): JSX.Element {
  return (
    <kbd
      data-testid="kbkey"
      className={cx(
        "inline-block font-mono text-xs bg-plate rounded-sm px-1.5 py-px text-ink-soft",
        // Uniform hair-strong border, just thicker (2px) at the bottom edge — the
        // raised-key 3D affordance, NOT a colored accent stripe (checklist G).
        "border border-hair-strong border-b-2",
        className,
      )}
    >
      {children}
    </kbd>
  );
}
