import type { ReactNode } from "react";
import { Button, Dot } from "../ui";

/**
 * Verdict — the loupe's tri-state screening verdict (SA-SCREENED).
 *
 * Three honest states: `dropped` (rejected), `kept` (explicitly accepted),
 * neither (unscreened). The old binary read "Kept" for a null status, which
 * made diligent screening of clean samples invisible — Keep is now a verb.
 *
 * Two toggles mirror the keyboard exactly:
 * - the drop toggle (`onToggle`, key X): Drop ↔ Restore (rejected/null)
 * - the keep toggle (`onToggleKeep`, key K): Keep ↔ Restore (accepted/null)
 * Last verb wins: Drop on a kept frame rejects directly; Keep on a dropped
 * frame accepts directly. No forced trip through unscreened.
 */
export interface VerdictProps {
  /** true = exposure dropped/rejected. */
  dropped: boolean;
  /** true = exposure explicitly kept/accepted. Ignored while `dropped`. */
  kept?: boolean;
  /** Optional hint override (default copy depends on the state). */
  hint?: ReactNode;
  /** Toggle drop/restore (the X verb). */
  onToggle?: () => void;
  /** Toggle keep/restore (the K verb). */
  onToggleKeep?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function Verdict({
  dropped,
  kept = false,
  hint,
  onToggle,
  onToggleKeep,
  className,
}: VerdictProps): JSX.Element {
  const isKept = !dropped && kept;
  const word = dropped ? "Dropped" : isKept ? "Kept" : "Unscreened";
  const defaultHint = dropped
    ? "Restore to clear the call."
    : isKept
      ? "Drop or restore to change the call."
      : "Keep or drop to screen this exposure.";
  const tone = dropped ? "accent" : isKept ? "success" : "neutral";

  return (
    <div
      className={`flex items-center gap-2.5 rounded-md border border-hair-strong bg-paper-sunk px-3 py-2.5${className ? ` ${className}` : ""}`}
    >
      <Dot tone={tone} size="sm" aria-hidden />

      <div className="flex-1">
        <div className="text-body font-bold text-ink">{word}</div>
        <div className="text-caption text-ink-soft">
          {hint ?? defaultHint}
        </div>
      </div>

      <Button
        variant="outline"
        {...(onToggleKeep ? { onClick: onToggleKeep } : {})}
      >
        {isKept ? "Restore" : "Keep"}
      </Button>
      <Button
        variant="outline"
        {...(onToggle ? { onClick: onToggle } : {})}
      >
        {dropped ? "Restore" : "Drop"}
      </Button>
    </div>
  );
}
