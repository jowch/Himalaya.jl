import { phaseColor } from "../phases";

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
    const color = phaseColor(phase);
    return (
      <span
        data-testid="phase-chip"
        className="text-data-strong rounded px-2 py-0.5 text-[11px] font-bold"
        style={{
          color,
          background: `color-mix(in oklab, ${color} 13%, transparent)`,
        }}
      >
        {phase}
      </span>
    );
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
