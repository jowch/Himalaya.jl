import type { HTMLAttributes } from "react";

export type DotTone = "accent" | "success" | "muted" | "neutral";
export type DotSize = "xs" | "sm";

interface DotProps extends Omit<HTMLAttributes<HTMLSpanElement>, "color"> {
  label?: string;
  tone?: DotTone;
  size?: DotSize;
}

const toneClass: Record<DotTone, string> = {
  accent: "bg-print-accent",
  success: "bg-success",
  neutral: "bg-hair-strong",
  muted: "border-[1.5px] border-hair-strong",
};

const sizeClass: Record<DotSize, string> = {
  xs: "h-1.5 w-1.5",
  sm: "h-2 w-2",
};

export function Dot({
  label,
  tone = "neutral",
  size = "sm",
  className = "",
  ...props
}: DotProps): JSX.Element {
  const decorative = props["aria-hidden"] === true || props["aria-hidden"] === "true" || label == null;
  return (
    <span
      className={`inline-block shrink-0 rounded-full ${sizeClass[size]} ${toneClass[tone]} ${className}`}
      data-tone={tone}
      {...(decorative ? {} : { role: "img", "aria-label": label })}
      {...props}
    />
  );
}
