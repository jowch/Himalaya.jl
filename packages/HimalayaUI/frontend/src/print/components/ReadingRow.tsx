import type { ReactNode } from "react";
import { Swatch, PhaseLabel } from "../ui";

export interface ReadingRowProps {
  phase: string;
  span: ReactNode; // e.g. "1:0 → 1:0.75"
  lattice: ReactNode; // e.g. "a 205 → 195 Å"
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

export function ReadingRow({ phase, span, lattice, className }: ReadingRowProps): JSX.Element {
  return (
    <div
      data-testid="reading-row"
      className={`flex flex-col gap-0.5${className ? ` ${className}` : ""}`}
    >
      <div className="flex items-center gap-[7px]">
        <Swatch phase={phase} shape="circle" />
        <PhaseLabel phase={phase} className="text-meta font-bold">
          {phase}
        </PhaseLabel>
        <span className="ml-auto text-data text-ink-soft">{span}</span>
      </div>
      <div className="text-data text-ink-faint ml-4">{lattice}</div>
    </div>
  );
}
