import { phaseColor } from "../../phases";

export function PhaseChip({ phase }: { phase: string }): JSX.Element {
  return (
    <span
      className="inline-flex items-center rounded-sm px-1.5 text-xs font-semibold leading-[1.6] text-plate"
      style={{ background: phaseColor(phase) }}
    >
      {phase}
    </span>
  );
}
