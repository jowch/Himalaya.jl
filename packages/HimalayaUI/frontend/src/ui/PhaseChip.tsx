import { Fragment, type CSSProperties, type HTMLAttributes } from "react";
import { phaseColor } from "../phases";
import { PhaseShortName } from "./PhaseName";

export type PhaseChipVariant = "tint" | "solid" | "text";
export type PhaseChipSize = "sm" | "md";


interface PhaseChipProps
  extends Omit<HTMLAttributes<HTMLSpanElement>, "color" | "children" | "title"> {
  phase?: string;
  formFactor?: boolean;
  coexistWith?: string[];
  variant?: PhaseChipVariant;
  size?: PhaseChipSize;
  className?: string;
}
const base = "inline-flex items-center rounded-sm border";
const sizeClass: Record<PhaseChipSize, string> = {
  sm: "text-xs font-semibold px-2 py-0.5",
  md: "text-base font-semibold px-1.5 py-0.5",
};
const textSize: Record<PhaseChipSize, string> = {
  sm: "text-xs font-semibold",
  md: "text-base font-semibold",
};
function variantStyle(variant: PhaseChipVariant, color: string): CSSProperties {
  if (variant === "text") return { color };
  if (variant === "solid") {
    return {
      color: "var(--color-paper)",
      background: color,
      borderColor: "transparent",
    };
  }
  return {
    color,
    background: `color-mix(in oklab, ${color} 13%, transparent)`,
    borderColor: `color-mix(in oklab, ${color} 35%, transparent)`,
  };
}
export function PhaseChip({
  phase,
  coexistWith,
  formFactor = false,
  variant = "tint",
  size = "sm",
  className = "",
  ...props
}: PhaseChipProps): JSX.Element {
  if (formFactor) {
    return (
      <span
        data-testid="phase-chip"
        data-form-factor="true"
        data-size={size}
        className={`${base} ${sizeClass[size]} ${className}`}
        style={{
          color: "var(--color-ink-soft)",
          background: "var(--color-paper-sunk)",
          borderColor: "var(--color-hair-strong)",
        }}
        {...props}
      >
        Form factor
      </span>
    );
  }
  const p = phase ?? "";
  const color = phaseColor(p);
  const coexist = coexistWith ?? [];
  const text =
    coexist.length > 0
      ? [p, ...coexist].map((q, i) => (
          <Fragment key={q}>
            {i > 0 && " + "}
            <PhaseShortName phase={q} />
          </Fragment>
        ))
      : p;
  return (
    <span
      data-testid="phase-chip"
      data-variant={variant}
      data-size={size}
      data-coexist={coexist.length > 0 ? "true" : undefined}
      data-coexist-count={coexist.length > 0 ? String(coexist.length + 1) : undefined}
      className={variant === "text" ? `inline-flex items-center ${textSize[size]} ${className}` : `${base} ${sizeClass[size]} ${className}`}
      style={variantStyle(variant, color)}
      {...props}
    >
      {text}
    </span>
  );
}