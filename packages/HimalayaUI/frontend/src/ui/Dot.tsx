import type { HTMLAttributes } from "react";

export type DotTone = "accent" | "success" | "warning" | "muted" | "neutral";
export type DotSize = "xs" | "sm" | "md";

interface DotProps extends Omit<HTMLAttributes<HTMLSpanElement>, "color" | "title"> {
  label?: string;
  tone?: DotTone;
  size?: DotSize;
  bordered?: boolean;
}
const toneClass: Record<DotTone, string> = {
  accent: "bg-print-accent",
  success: "bg-success",
  warning: "bg-warning",
  neutral: "bg-hair-strong",
  muted: "border-[1.5px] border-hair-strong",
};
const sizeClass: Record<DotSize, string> = {
  xs: "h-1.5 w-1.5",
  sm: "h-2 w-2",
  md: "h-[9px] w-[9px]",
};
export function Dot({
  label,
  tone = "neutral",
  size = "sm",
  bordered = false,
  className = "",
  ...props
}: DotProps): JSX.Element {
  const decorative = props["aria-hidden"] === true || props["aria-hidden"] === "true" || label == null;
  const borderClass = bordered ? " border-[1.5px] border-plate" : "";
  return (
    <span
      className={`inline-block shrink-0 rounded-full ${sizeClass[size]} ${toneClass[tone]}${borderClass} ${className}`}
      data-tone={tone}
      data-size={size}
      {...(bordered ? { "data-bordered": "true" } : {})}
      {...(decorative ? {} : { role: "img", "aria-label": label })}
      {...props}
    />
  );
}