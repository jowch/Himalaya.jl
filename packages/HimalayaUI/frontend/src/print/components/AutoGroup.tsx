import type { ReactNode } from "react";

export interface AutoGroupProps {
  /** Body copy — ReactNode so the caller can embed <strong> emphasis. */
  children: ReactNode;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

function Star(): JSX.Element {
  return (
    <svg
      data-role="autogroup-star"
      className="w-[15px] h-[15px] shrink-0"
      viewBox="0 0 16 16"
      fill="none"
      aria-hidden="true"
    >
      <path
        d="M8 1.4l1.7 4.2 4.5.3-3.5 2.9 1.2 4.4L8 10.9 4.1 13.2l1.2-4.4L1.8 5.9l4.5-.3z"
        fill="var(--color-print-accent)"
      />
    </svg>
  );
}

export function AutoGroup({ children, className }: AutoGroupProps): JSX.Element {
  return (
    <div
      data-testid="auto-group"
      data-variant="summary"
      className={`rounded border border-hair p-3 bg-paper-sunk${className ? ` ${className}` : ""}`}
    >
      <div className="flex gap-2 items-start">
        <Star />
        <div className="text-body text-ink-soft">{children}</div>
      </div>
    </div>
  );
}
