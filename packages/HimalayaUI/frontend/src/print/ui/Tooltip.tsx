import { cloneElement, useId, useState } from "react";
import type { ReactElement } from "react";
import { cx } from "../../lib/cx";

interface TooltipProps {
  /** The caption text shown in the popover. */
  label: string;
  /** A single focusable/hoverable trigger element. */
  children: ReactElement;
  /** Which side of the trigger the tip sits on. Default `"top"`. */
  side?: "top" | "bottom";
  /** Allow the caption to wrap (capped width) instead of a single nowrap strip.
   *  For multi-clause guidance that would otherwise overflow the viewport. */
  multiline?: boolean;
  /** Placement-only consumer override on the tip span (closed-look contract). */
  className?: string;
}


/**
 * Caption popover that shows on hover AND focus around a single trigger child.
 *
 * The Tooltip is the intentional DARK exception in the warm-paper system
 * (DESIGN.md §2 frame-edge / §4 detector dark window): a dark `frame-edge`
 * window with light `frame-tag` caption text, set via inline token vars so the
 * design guard's dark-exception stays explicit. A11y: `role="tooltip"` plus
 * `aria-describedby` wired from the trigger to the tip id (React `useId()`,
 * no `Math.random()`). No entrance bounce — opacity only; reduced-motion is
 * handled globally.
 */
export function Tooltip({
  label,
  children,
  side = "top",
  multiline = false,
  className = "",
}: TooltipProps): JSX.Element {
  const id = useId();
  const [open, setOpen] = useState(false);

  // Only add aria-describedby when open, to avoid passing explicit `undefined`
  // under exactOptionalPropertyTypes (conditional spread).
  const trigger = cloneElement(
    children,
    open ? { "aria-describedby": id } : {},
  );

  return (
    <span
      className="relative inline-flex"
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => setOpen(false)}
      onFocus={() => setOpen(true)}
      onBlur={() => setOpen(false)}
    >
      {trigger}
      {open && (
        <span
          role="tooltip"
          id={id}
          data-testid="tooltip"
          data-side={side}
          className={cx(
            "pointer-events-none absolute z-30 left-1/2 -translate-x-1/2 rounded-sm px-2.5 py-1 text-xs font-medium shadow-lg",
            multiline ? "w-max max-w-[15rem] whitespace-normal text-left leading-snug" : "whitespace-nowrap",
            side === "top" ? "bottom-full mb-1.5" : "top-full mt-1.5",
            className,
          )}
          style={{
            background: "var(--color-frame-edge)",
            color: "var(--color-frame-tag)",
            // A hair of inset definition so the dark window reads against dark
            // detector panels too (the warm-paper system's only dark exception).
            boxShadow: "0 6px 16px -4px rgba(0,0,0,0.35), inset 0 0 0 1px rgba(255,255,255,0.06)",
          }}
        >
          {label}
        </span>
      )}
    </span>
  );
}
