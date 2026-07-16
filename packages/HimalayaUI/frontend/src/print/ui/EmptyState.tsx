import type { ReactNode } from "react";
import { cx } from "../../lib/cx";


interface EmptyStateProps {
  title: string;
  body?: ReactNode;
  /** A way forward — the consumer composes a ui control (e.g. a Button wired
   *  to a refetch or navigation); the primitive owns the gap below the body. */
  action?: ReactNode;
  /** Heading level for the title. Default "h2" — EmptyState usually sits inside
   *  a page that already owns the h1. On a route where EmptyState IS the page
   *  (a not-found view with no other heading), pass "h1" so the document is not
   *  left with no top-level heading / a level-skip (FO-NFHEAD, WCAG 1.3.1). The
   *  visual treatment is identical at either level — level is semantics, not
   *  size. */
  as?: "h1" | "h2";
  /** Placement-only className (spacing, max-width, etc). Appearance is fixed. */
  className?: string;
}

export function EmptyState({ title, body, action, as = "h2", className }: EmptyStateProps): JSX.Element {
  const Heading = as;
  return (
    <div data-testid="empty-state" className={cx("text-center py-16 px-5", className)}>
      {/* 21px serif in the mockup; reuse the .text-headline role (19px, -2px) rather
          than mint a one-off token (YAGNI; the audit warns against off-scale sizes). */}
      <Heading className="text-headline text-ink-soft mb-1.5">{title}</Heading>
      {body && <div className="text-base text-ink-soft">{body}</div>}
      {action && <div className="mt-4 flex justify-center">{action}</div>}
    </div>
  );
}
