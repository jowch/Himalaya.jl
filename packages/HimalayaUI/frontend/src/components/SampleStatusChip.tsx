import { PhaseChip } from "./ui/PhaseChip";

/**
 * The contact-sheet status cell (M-6). A sample with a resolved phase shows a
 * monospace phase chip tinted in the R0b phase colour (consumed via
 * `phaseColor`, never re-defined here); a sample without one shows the
 * hollow-dot "Not indexed" affordance from sample-table.html's `.unset`.
 *
 * Colour is never the sole signal — the chip always carries the phase name
 * (DESIGN.md colour-blind note).
 */
export function SampleStatusChip({
  phase,
}: {
  phase: string | null | undefined;
}): JSX.Element {
  if (phase) {
    return <PhaseChip phase={phase} variant="tint" size="sm" />;
  }
  return (
    <span className="flex items-center text-xs text-ink-faint">
      <span
        data-testid="status-dot"
        className="mr-1.5 inline-block h-1.5 w-1.5 rounded-full border-[1.5px] border-hair-strong"
      />
      Not indexed
    </span>
  );
}
