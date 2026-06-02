import { phaseColor } from "../../phases";

export type SwatchShape = "square" | "circle";

interface SwatchProps {
  /** Phase whose color fills the swatch (via phaseColor). */
  phase: string;
  /** "square" (rounded color chip, default) | "circle" (reading-row dot). */
  shape?: SwatchShape;
  /** PLACEMENT-ONLY. */
  className?: string;
}

const shapeClass: Record<SwatchShape, string> = {
  square: "rounded-sm",
  circle: "rounded-full",
};

export function Swatch({ phase, shape = "square", className = "" }: SwatchProps): JSX.Element {
  return (
    <span
      data-swatch
      data-phase={phase}
      data-shape={shape}
      aria-hidden="true"
      className={`inline-block shrink-0 h-[9px] w-[9px] ${shapeClass[shape]} ${className}`}
      style={{ background: phaseColor(phase) }}
    />
  );
}
// later: gradient/coexistence mode (series-plot member)
