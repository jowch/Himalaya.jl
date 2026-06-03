import { phaseColor } from "../../phases";

export type SwatchShape = "square" | "circle";
export type SwatchSize = "sm" | "md";

interface SwatchProps {
  /** Phase whose color fills the swatch (via phaseColor). Ignored when `empty`. */
  phase: string;
  /** "square" (rounded color chip, default) | "circle" (reading-row dot). */
  shape?: SwatchShape;
  /** "sm" 9px (default) | "md" 11px (series-plot member). */
  size?: SwatchSize;
  /** Second phase → a 135° gradient blending both (coexistence). */
  coexistWith?: string;
  /** Form factor: transparent fill + dashed hairline, no phase color. */
  empty?: boolean;
  /** PLACEMENT-ONLY. */
  className?: string;
}

const shapeClass: Record<SwatchShape, string> = {
  square: "rounded-sm",
  circle: "rounded-full",
};
const sizeClass: Record<SwatchSize, string> = {
  sm: "h-[9px] w-[9px]",
  md: "h-[11px] w-[11px]",
};

export function Swatch({
  phase,
  shape = "square",
  size = "sm",
  coexistWith,
  empty = false,
  className = "",
}: SwatchProps): JSX.Element {
  const background = empty
    ? "transparent"
    : coexistWith
      ? `linear-gradient(135deg, ${phaseColor(phase)} 48%, ${phaseColor(coexistWith)} 52%)`
      : phaseColor(phase);
  return (
    <span
      data-swatch
      data-phase={phase}
      data-shape={shape}
      {...(coexistWith ? { "data-coexist": coexistWith } : {})}
      {...(empty ? { "data-empty": "true" } : {})}
      aria-hidden="true"
      className={`inline-block shrink-0 ${sizeClass[size]} ${shapeClass[shape]} ${empty ? "border-[1.4px] border-dashed border-hair-strong" : ""} ${className}`.trim()}
      style={{ background }}
    />
  );
}
