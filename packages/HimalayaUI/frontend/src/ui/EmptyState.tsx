import type { ReactNode } from "react";
import { cx } from "../lib/cx";

interface EmptyStateProps {
  title: string;
  body?: ReactNode;
  action?: ReactNode;
  as?: "h1" | "h2";
  className?: string;
}
export function EmptyState({ title, body, action, as = "h2", className }: EmptyStateProps): JSX.Element {
  const Heading = as;
  return (
    <div data-testid="empty-state" className={cx("text-center py-16 px-5", className)}>
      <Heading className="text-headline text-ink-soft mb-1.5">{title}</Heading>
      {body && <div className="text-base text-ink-soft">{body}</div>}
      {action && <div className="mt-4 flex justify-center">{action}</div>}
    </div>
  );
}