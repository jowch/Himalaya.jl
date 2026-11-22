import type { CSSProperties, HTMLAttributes } from "react";
import { phaseColor } from "../phases";

export type PhaseChipVariant = "tint" | "solid";
export type PhaseChipSize = "sm" | "md";

const PHASE_SHORT: Record<string, string> = {
  Pn3m: "Pn3m", Im3m: "Im3m", Ia3d: "Ia3d", Fm3m: "Fm3m",
  Fd3m: "Fd3m", Hexagonal: "Hex", Lamellar: "Lam", Square: "Sq",
};
function short(p: string): string { return PHASE_SHORT[p] ?? p; }

interface PhaseChipProps
  extends Omit<HTMLAttributes<HTMLSpanElement>, "color" | "children"> {
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
function variantStyle(variant: PhaseChipVariant, color: string): CSSProperties {
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
      ? [p, ...coexist].map(short).join(" + ")
      : p;
  return (
    <span
      data-testid="phase-chip"
      data-variant={variant}
      data-size={size}
      data-coexist={coexist.length > 0 ? "true" : undefined}
      data-coexist-count={coexist.length > 0 ? String(coexist.length + 1) : undefined}
      className={`${base} ${sizeClass[size]} ${className}`}
      style={variantStyle(variant, color)}
      {...props}
    >
      {text}
    </span>
  );
}