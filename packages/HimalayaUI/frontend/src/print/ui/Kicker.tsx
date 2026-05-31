import type { HTMLAttributes } from "react";

export type KickerTone = "accent" | "faint";

interface KickerProps extends HTMLAttributes<HTMLElement> {
  /** terracotta page/section eyebrow vs. dim metric/field label. Default: "faint". */
  tone?: KickerTone;
  /**
   * Element to render. Section headings that introduce a region should be a
   * heading ("h2"/"h3") for the a11y tree; inline/aside labels stay "div"/"span".
   * Default: "div" (matches the prevailing inline usage).
   */
  as?: "div" | "span" | "h2" | "h3";
}

// tone selects color + per-tone tracking; base geometry (700/uppercase/size) is .text-kicker.
const toneClass: Record<KickerTone, string> = {
  accent: "text-kicker-accent",
  faint: "text-kicker-faint",
};

export function Kicker({
  tone = "faint",
  as: Tag = "div",
  className = "",
  children,
  ...props
}: KickerProps): JSX.Element {
  return (
    <Tag
      data-tone={tone}
      className={`text-kicker ${toneClass[tone]} ${className}`}
      {...props}
    >
      {children}
    </Tag>
  );
}
