import type { ReactNode } from "react";
import { phaseColor } from "../../phases";

interface PhaseLabelProps {
  /** Phase whose color tints the text (via phaseColor). */
  phase: string;
  as?: "span" | "div";
  /** PLACEMENT-ONLY (font/size/weight live here — e.g. text-display, or text-data-strong font-bold). */
  className?: string;
  children: ReactNode;
}

export function PhaseLabel({ phase, as: Tag = "span", className = "", children }: PhaseLabelProps): JSX.Element {
  return (
    <Tag data-phase-label data-phase={phase} className={className} style={{ color: phaseColor(phase) }}>
      {children}
    </Tag>
  );
}
