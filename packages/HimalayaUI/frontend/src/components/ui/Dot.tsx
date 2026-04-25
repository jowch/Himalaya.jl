import type { HTMLAttributes } from "react";

interface DotProps extends HTMLAttributes<HTMLSpanElement> {
  label: string;
}

export function Dot({ label, className = "", ...props }: DotProps): JSX.Element {
  return (
    <span
      className={`inline-block w-2 h-2 rounded-full flex-shrink-0 ${className}`}
      aria-label={label}
      role="img"
      {...props}
    />
  );
}
