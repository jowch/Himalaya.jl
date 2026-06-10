import type { ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface EmptyStateProps {
  title: string;
  body?: ReactNode;
  /** A way forward — the consumer composes a ui control (e.g. a Button wired
   *  to a refetch or navigation); the primitive owns the gap below the body. */
  action?: ReactNode;
  /** Placement-only className (spacing, max-width, etc). Appearance is fixed. */
  className?: string;
}

export function EmptyState({ title, body, action, className }: EmptyStateProps): JSX.Element {
  return (
    <div data-testid="empty-state" className={cx("text-center py-16 px-5", className)}>
      {/* 21px serif in the mockup; reuse the .text-headline role (19px, -2px) rather
          than mint a one-off token (YAGNI; the audit warns against off-scale sizes). */}
      <h2 className="text-headline text-ink-soft mb-1.5">{title}</h2>
      {body && <div className="text-base text-ink-soft">{body}</div>}
      {action && <div className="mt-4 flex justify-center">{action}</div>}
    </div>
  );
}
