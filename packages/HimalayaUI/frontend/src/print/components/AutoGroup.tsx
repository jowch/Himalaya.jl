import type { ReactNode } from "react";

export interface AutoGroupAction {
  label: string;
  onClick?: () => void;
  muted?: boolean;
}

export interface AutoGroupProps {
  /** "summary" (recessed, scoping) | "compose" (plate bg, builder). Default "summary". */
  variant?: "summary" | "compose";
  /** Optional bold title (the compose variant's "Auto-grouped"). */
  title?: string;
  /** Body copy — ReactNode so the caller can embed <strong> emphasis. */
  children: ReactNode;
  /** Optional link-style actions (compose variant). */
  actions?: AutoGroupAction[];
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

export function AutoGroup({
  variant = "summary",
  title,
  children,
  actions,
  className,
}: AutoGroupProps): JSX.Element {
  const bg = variant === "compose" ? "bg-plate" : "bg-paper-sunk";
  const body = <div className="text-body text-ink-soft">{children}</div>;
  return (
    <div
      data-testid="auto-group"
      data-variant={variant}
      className={`rounded border border-hair p-3 ${bg}${className ? ` ${className}` : ""}`}
    >
      <div className="flex gap-2 items-start">
        <Star />
        {title ? (
          <span className="text-meta text-ink font-bold">{title}</span>
        ) : (
          body
        )}
      </div>

      {title && <div className="mt-2">{body}</div>}

      {actions?.length ? (
        <div className="flex gap-3 mt-2">
          {actions.map((action) => (
            <button
              key={action.label}
              type="button"
              onClick={action.onClick}
              className={
                action.muted
                  ? "text-sm font-semibold text-ink-faint hover:underline"
                  : "text-sm font-semibold text-print-accent hover:underline"
              }
            >
              {action.label}
            </button>
          ))}
        </div>
      ) : null}
    </div>
  );
}
