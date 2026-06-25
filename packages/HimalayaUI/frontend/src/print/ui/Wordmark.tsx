import type { ReactNode } from "react";
import { cx } from "../../lib/cx";

interface WordmarkProps {
  children: ReactNode;
  /** Optional faint tail after a middle dot, e.g. "SAXS" -> "HIMALAYA · SAXS". */
  tail?: ReactNode;
  className?: string;
}


/** The brand wordmark (mockup `.wordmark`): SANS, uppercase, weight 700,
 *  wide tracking — NOT a serif title (serif is reserved for plate titles).
 *  `tracking-[0.16em]` is the brand mark's specified tracking (an intentional
 *  in-print/ui arbitrary value for this one signature element); requires the
 *  Plus Jakarta Sans 700 weight (imported in styles.css). */
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
