import type { ReactNode } from "react";
import { Button, Dot } from "../ui";

/**
 * Verdict — the loupe's binary screening verdict.
 *
 * A frame is either kept (status null, the default) or dropped (status
 * rejected). There is one verb: Drop — a TOGGLE. Dropping a kept frame culls
 * it; pressing Drop again brings it back. There is no Keep / Restore action;
 * "Kept" is just the un-culled state. State shows via the word + dot; the Drop
 * button reflects the toggle with `aria-pressed`.
 */
export interface VerdictProps {
  /** true = exposure dropped (status rejected); false = kept (status null). */
  dropped: boolean;
  /** Optional hint override (default copy depends on the state). */
  hint?: ReactNode;
  /** Drop toggle (key X): rejected ↔ null. */
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
  const word = dropped ? "Dropped" : "Kept";
  const defaultHint = dropped
    ? "Drop again to bring this frame back." // LO-TERM: "frame", not "exposure"
    : "Drop to cull this frame.";
  const tone = dropped ? "accent" : "success";

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
        variant={dropped ? "accent" : "outline"}
        aria-pressed={dropped}
        {...(onToggle ? { onClick: onToggle } : {})}
      >
        Drop
      </Button>
    </div>
  );
}
