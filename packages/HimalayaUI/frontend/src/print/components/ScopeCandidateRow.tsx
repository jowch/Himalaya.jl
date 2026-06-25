import type { ReactNode } from "react";
import { Sparkline } from "../plot/Sparkline";
import { Button } from "../ui";
import { cx } from "../../lib/cx";


export interface ScopeCandidateRowProps {
  name: string;
  /** Mono `smp_{id}` identity, mirroring the member-row (ScopeSampleRow)
   *  idiom — display names duplicate in real corpora. */
  sampleId?: string;
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
  sampleId,
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
        <div className="flex items-center gap-2 min-w-0">
          {/* SC-CANDNAME-DIM: name at text-ink (primary identifier); the row's
              secondary status comes from the dashed border + dimmed sparkline. */}
          <span className="text-meta font-semibold text-ink">{name}</span>
          {sampleId != null && (
            <span className="text-caption font-mono text-ink-soft">{sampleId}</span>
          )}
        </div>
        <div className="text-caption text-ink-soft">{why}</div>
      </div>
      <Button variant="outline" {...(onAdd ? { onClick: onAdd } : {})}>
        + Add to series
      </Button>
    </div>
  );
}
