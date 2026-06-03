import type { ReactNode } from "react";
import { Card } from "../ui";
import { ReadingRow } from "./ReadingRow";

export interface ReadingDatum {
  phase: string;
  span: ReactNode;
  lattice: ReactNode;
}

export interface ReadingPanelProps {
  readings: ReadingDatum[];
  /** e.g. "coexistence at 1:0.5". */
  coexistenceNote?: ReactNode;
  /** e.g. "form factor only at 1:1.5". */
  formFactorNote?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function ReadingPanel({ readings, coexistenceNote, formFactorNote, className }: ReadingPanelProps): JSX.Element {
  return (
    <Card data-testid="reading-panel" className={`px-[13px] py-[11px] flex flex-col gap-[9px]${className ? ` ${className}` : ""}`}>
      {readings.map((r, i) => (
        <ReadingRow key={i} phase={r.phase} span={r.span} lattice={r.lattice} />
      ))}
      {coexistenceNote != null && (
        <div data-testid="reading-coex" className="text-data text-ink-soft mt-0.5 pt-2 border-t border-hair">{coexistenceNote}</div>
      )}
      {formFactorNote != null && (
        <div data-testid="reading-ff" className="text-data text-ink-soft mt-0.5 pt-2 border-t border-hair">{formFactorNote}</div>
      )}
    </Card>
  );
}
