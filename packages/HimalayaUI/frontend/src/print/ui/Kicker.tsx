import type { HTMLAttributes } from "react";

export type KickerTone = "accent" | "faint" | "soft";

interface KickerProps extends HTMLAttributes<HTMLElement> {
  /**
   * "accent": terracotta page/section eyebrow.
   * "faint": dim decorative/redundant label — reserved for plain paper
   * (ink-faint clears only the 3:1 large/non-text line there, and fails
   * outright on sunk surfaces).
   * "soft": informational label, and ANY kicker on a sunk surface
   * (ink-soft, AA at all sizes — WCAG 1.4.3).
   * Default: "faint".
   */
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
  soft: "text-kicker-soft",
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
