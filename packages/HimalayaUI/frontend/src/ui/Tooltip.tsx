import { cloneElement, useId, useState } from "react";
import type { ReactElement } from "react";
import { cx } from "../lib/cx";

interface TooltipProps {
  label: string;
  children: ReactElement;
  side?: "top" | "bottom";
  multiline?: boolean;
  className?: string;
}

export function Tooltip({
  label,
  children,
  side = "top",
  multiline = false,
  className = "",
}: TooltipProps): JSX.Element {
  const id = useId();
  const [open, setOpen] = useState(false);
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
            boxShadow: "0 6px 16px -4px rgba(0,0,0,0.35), inset 0 0 0 1px rgba(255,255,255,0.06)",
          }}
        >
          {label}
        </span>
      )}
    </span>
  );
}