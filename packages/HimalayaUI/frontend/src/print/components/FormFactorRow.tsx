import { Card, CheckCircle } from "../ui";

export interface FormFactorRowProps {
  /** Whether the sample is currently declared form factor. */
  selected: boolean;
  /** Toggle the form-factor declaration on/off. */
  onToggle: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * FormFactorRow — the "this sample has no Bragg peaks to index" declaration.
 *
 * Pinned at the foot of the candidate list, it mirrors {@link CandidateRow}'s
 * shape (a selectable Card with a CheckCircle mark) so declaring form factor
 * feels like the same one-check interaction as picking a phase — but it carries
 * no PhaseChip or score, because it is not a phase. It is mutually exclusive
 * with the phase candidates: the backend clears any members when the assignment
 * goes to `form_factor`, and checking a phase flips the assignment back to
 * `indexed`. The page only renders this row while nothing is in the call, so a
 * declaration can never silently discard a confirmed phase.
 */
export function FormFactorRow({
  selected,
  onToggle,
  className,
}: FormFactorRowProps): JSX.Element {
  return (
    <Card
      as="button"
      selected={selected}
      aria-pressed={selected}
      aria-label={`Form factor, no Bragg peaks${selected ? ", declared" : ""}`}
      data-testid="form-factor-row"
      className={`w-full flex items-center gap-3 px-3 py-2.5 text-left${className ? ` ${className}` : ""}`}
      onClick={onToggle}
    >
      <CheckCircle checked={selected} />
      <div className="flex-1 min-w-0">
        <div className="text-body text-ink">Form factor</div>
        <div className="text-caption text-ink-soft mt-0.5">No Bragg peaks to index</div>
      </div>
    </Card>
  );
}
