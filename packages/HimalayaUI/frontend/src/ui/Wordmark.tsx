import type { ReactNode } from "react";
import { cx } from "../lib/cx";

interface WordmarkProps {
  children: ReactNode;
  tail?: ReactNode;
  className?: string;
}

export function Wordmark({ children, tail, className = "" }: WordmarkProps): JSX.Element {
  return (
    <span
      data-testid="wordmark"
      className={cx("text-lg font-semibold text-ink", className)}
    >
      {children}
      {tail && (
        <span data-role="tail" className="text-ink-soft font-semibold"> · {tail}</span>
      )}
    </span>
  );
}
