import type { ReactNode } from "react";
import { Button, Dot } from "../ui";

export interface VerdictProps {
  /** true = exposure dropped/rejected; false = kept. */
  dropped: boolean;
  /** Optional hint override (default copy depends on dropped). */
  hint?: ReactNode;
  /** Toggle drop/restore. */
  onToggle?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function Verdict({
  dropped,
  hint,
  onToggle,
  className,
}: VerdictProps): JSX.Element {
  const defaultHint = dropped
    ? "Restore to keep this exposure."
    : "Everything is kept until you drop it.";

  return (
    <div
      className={`flex items-center gap-2.5 rounded-md border border-hair-strong bg-paper-sunk px-3 py-2.5${className ? ` ${className}` : ""}`}
    >
      <Dot tone={dropped ? "accent" : "success"} size="sm" aria-hidden />

      <div className="flex-1">
        <div className="text-body font-bold text-ink">
          {dropped ? "Dropped" : "Kept"}
        </div>
        <div className="text-caption text-ink-faint">
          {hint ?? defaultHint}
        </div>
      </div>

      <Button
        variant="outline"
        {...(onToggle ? { onClick: onToggle } : {})}
      >
        {dropped ? "Restore" : "Drop"}
      </Button>
    </div>
  );
}
