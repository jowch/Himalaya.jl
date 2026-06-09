import type { ReactNode } from "react";
import { Sparkline } from "../plot/Sparkline";
import { Button } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ScopeCandidateRowProps {
  name: string;
  /** One-line rationale; ReactNode so the caller can embed an accent
   *  `<strong>` for the `.warn` emphasis. */
  why: ReactNode;
  trace: { q: number[]; I: number[] };
  phase?: string | null;
  onAdd?: () => void;
  /** PLACEMENT ONLY. */
  className?: string;
}

/** A scoping sample-candidate row (the mockup `.cand-row`): a dashed,
 *  hairline-bordered, rounded row with a dimmed sparkline, the candidate
 *  name + a one-line "why", and a trailing "+ Add to series" button. */
export function ScopeCandidateRow({
  name,
  why,
  trace,
  phase,
  onAdd,
  className,
}: ScopeCandidateRowProps): JSX.Element {
  return (
    <div
      data-testid="scope-candidate-row"
      className={cx(
        "flex items-center gap-3 px-2 py-2.5 border border-dashed border-hair-strong rounded",
        className,
      )}
    >
      <Sparkline
        trace={trace}
        {...(phase != null ? { phase } : {})}
        className="opacity-70"
      />
      <div className="flex-1 min-w-0">
        <div className="text-meta font-semibold text-ink-soft">{name}</div>
        <div className="text-caption text-ink-soft">{why}</div>
      </div>
      <Button variant="outline" {...(onAdd ? { onClick: onAdd } : {})}>
        + Add to series
      </Button>
    </div>
  );
}
