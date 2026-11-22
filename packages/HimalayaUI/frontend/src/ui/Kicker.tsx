import { forwardRef } from "react";
import type { HTMLAttributes } from "react";

export type KickerTone = "accent" | "faint" | "soft" | "warning";

interface KickerProps extends HTMLAttributes<HTMLElement> {
  tone?: KickerTone;
  as?: "div" | "span" | "h2" | "h3";
}
const toneClass: Record<KickerTone, string> = {
  accent: "text-kicker-accent",
  faint: "text-kicker-faint",
  soft: "text-kicker-soft",
  warning: "text-warning",
};
export const Kicker = forwardRef<HTMLElement, KickerProps>(function Kicker(
  { tone = "faint", as: Tag = "div", className = "", children, ...props },
  ref,
): JSX.Element {
  return (
    <Tag
      ref={ref as React.Ref<never>}
      data-tone={tone}
      className={`text-kicker ${toneClass[tone]} ${className}`}
      {...props}
    >
      {children}
    </Tag>
  );
});