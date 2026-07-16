import type { ReactNode } from "react";
import { GripHandle, Swatch, PhaseChip } from "../ui";

export interface MemberRowProps {
  name: string;
  sub?: ReactNode; // sample id / dose / variable value
  phase: string; // dominant phase (drives swatch + chip)
  coexistWith?: string[]; // optional coexisting phases for the chip
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

// ONE reorderable trace-list row in the series-builder sidebar (mockup `.trow`):
// grip + square color swatch + name/dose + phase chip. The root carries `group`
// so the grip's color reveal fires on row hover. Resting→hover cue is
// `hover:bg-plate` + a hairline border, matching the mockup. The border is
// reserved at rest as `border-transparent` so hover only changes its COLOR —
// adding a 1px border on hover would grow the box by 2px and shift the row.
export function MemberRow({
  name,
  sub,
  phase,
  coexistWith,
  className,
}: MemberRowProps): JSX.Element {
  return (
    <div
      data-testid="member-row"
      className={`group flex items-center gap-2 px-2 py-1.5 rounded border border-transparent hover:bg-plate hover:border-hair${className ? ` ${className}` : ""}`}
    >
      <GripHandle />
      <Swatch phase={phase} shape="square" />
      <div className="flex-1 min-w-0">
        <div className="text-meta font-semibold text-ink truncate">{name}</div>
        {sub != null && (
          <div className="text-data text-ink-soft truncate">{sub}</div>
        )}
      </div>
      <PhaseChip
        phase={phase}
        variant="tint"
        size="sm"
        className="shrink-0"
        {...(coexistWith ? { coexistWith } : {})}
      />
    </div>
  );
}
