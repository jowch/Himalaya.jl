import type { ReactNode } from "react";
import { cx } from "../../lib/cx";

interface WordmarkProps {
  children: ReactNode;
  /** Optional faint tail after a middle dot, e.g. "SAXS" -> "HIMALAYA · SAXS". */
  tail?: ReactNode;
  className?: string;
}

export function Wordmark({ children, tail, className = "" }: WordmarkProps): JSX.Element {
  return (
    <span
      data-testid="wordmark"
      className={cx("font-sans font-bold uppercase tracking-[0.16em] text-ink text-xs", className)}
    >
      {children}
      {tail && (
        <span data-role="tail" className="text-ink-soft font-semibold"> · {tail}</span>
      )}
    </span>
  );
}
