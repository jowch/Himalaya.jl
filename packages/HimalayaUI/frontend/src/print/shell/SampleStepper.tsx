import { IconButton, Kicker } from "../ui";

export interface SampleStepperProps {
  /** Active sample display name. */
  name: string;
  /** 0-based position of the active sample within its list. */
  index: number;
  /** Total samples in the list. */
  total: number;
  /** Step to the previous sample. Omit → the ‹ button is disabled. */
  onPrev?: () => void;
  /** Step to the next sample. Omit → the › button is disabled. */
  onNext?: () => void;
  /** Keyboard key that also steps prev (e.g. "["); shown as a tooltip + aria hint. */
  prevKey: string;
  /** Keyboard key that also steps next (e.g. "]"). */
  nextKey: string;
}

/**
 * The shared per-sample stepper for the TopBar's right slot: `‹ name / sample N
 * of M ›`. Used by BOTH the Focus workspace and the Loupe (same component, same
 * location) so inter-sample nav reads identically; each surface wires its own
 * prev/next (Focus steps experiment-siblings → /sample/:id, Loupe steps the
 * contact-sheet order → the loupe). The ‹ › buttons surface their `[`/`]`
 * keyboard equivalents via a tooltip (title) + aria-keyshortcuts.
 */
export function SampleStepper({
  name,
  index,
  total,
  onPrev,
  onNext,
  prevKey,
  nextKey,
}: SampleStepperProps): JSX.Element {
  return (
    <div data-testid="sample-stepper" className="flex items-center gap-2 text-ink">
      <IconButton
        label="Previous sample"
        title={`Previous sample (${prevKey})`}
        aria-keyshortcuts={prevKey}
        tone="ghost"
        disabled={onPrev === undefined}
        {...(onPrev ? { onClick: onPrev } : {})}
        data-testid="sample-stepper-prev"
      >
        {"‹"}
      </IconButton>
      <span className="flex flex-col items-end leading-tight">
        <span className="text-xs font-semibold text-ink">{name}</span>
        {/* F6: "sample N of M" is small INFORMATIONAL text → ink-soft (AA-normal);
            faint would fail AA at 11.5px. */}
        <Kicker as="span" tone="soft">
          sample {index + 1} of {total}
        </Kicker>
      </span>
      <IconButton
        label="Next sample"
        title={`Next sample (${nextKey})`}
        aria-keyshortcuts={nextKey}
        tone="ghost"
        disabled={onNext === undefined}
        {...(onNext ? { onClick: onNext } : {})}
        data-testid="sample-stepper-next"
      >
        {"›"}
      </IconButton>
    </div>
  );
}
