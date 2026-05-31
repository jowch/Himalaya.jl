import type { CSSProperties, HTMLAttributes } from "react";
import { phaseColor } from "../../phases";

export type PhaseChipVariant = "tint" | "solid";
export type PhaseChipSize = "sm" | "md";

const PHASE_SHORT: Record<string, string> = {
  Pn3m: "Pn3m", Im3m: "Im3m", Ia3d: "Ia3d", Fm3m: "Fm3m",
  Fd3m: "Fd3m", Hexagonal: "Hex", Lamellar: "Lam", Square: "Sq",
};
function short(p: string): string { return PHASE_SHORT[p] ?? p; }

interface PhaseChipProps
  extends Omit<HTMLAttributes<HTMLSpanElement>, "color" | "children"> {
  /** Phase name, e.g. "Pn3m". Rendered as the chip's text (the always-on
   *  second channel) AND drives the hue via phaseColor() internally. */
  phase: string;
  /** Optional second phase for a coexistence chip. When set, the chip reads
   *  "<short(phase)> + <short(coexistWith)>" using the DOMINANT (`phase`) color
   *  — a single tinted chip, never a split. (mockup builder coexistence chip) */
  coexistWith?: string;
  /** "tint": phase-tinted fill + phase-colored text (default, the M-6 look).
   *  "solid": phase-colored fill + paper text (emphasis / dense rows). */
  variant?: PhaseChipVariant;
  /** "sm" = 11px mono (contact sheet, inline rows). "md" = 13px (candidate
   *  rows, builder). Default "sm". */
  size?: PhaseChipSize;
  /** PLACEMENT ONLY — margin / width / grid-or-flex position. Appearance
   *  utilities are banned by the lint:design guard. */
  className?: string;
}

const base = "inline-flex items-center rounded-sm border";

const sizeClass: Record<PhaseChipSize, string> = {
  sm: "font-mono text-[11px] font-bold px-2 py-0.5",
  md: "font-mono text-[13px] font-bold px-1.5 py-0.5",
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
  variant = "tint",
  size = "sm",
  className = "",
  ...props
}: PhaseChipProps): JSX.Element {
  const color = phaseColor(phase);
  const text = coexistWith ? `${short(phase)} + ${short(coexistWith)}` : phase;
  return (
    <span
      data-testid="phase-chip"
      data-variant={variant}
      data-size={size}
      data-coexist={coexistWith ? "true" : undefined}
      className={`${base} ${sizeClass[size]} ${className}`}
      style={variantStyle(variant, color)}
      {...props}
    >
      {text}
    </span>
  );
}
