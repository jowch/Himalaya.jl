import { cx } from "../../lib/cx";
interface BonnetBadgeProps {
  /** PLACEMENT ONLY — margin / inline position. Appearance utilities are
   *  banned by the lint:design guard; the look lives in this file. */
  className?: string;
}


/** A small terracotta-tinted pill flagging the Gauss–Bonnet ratio signal on a
 *  candidate. Decorative status flag — non-interactive (no C/D interaction
 *  state or touch target; it carries no action).
 *
 *  Second-channel (checklist A): the flag is accent-tinted, but its meaning
 *  never rests on the tint alone — the literal WORD "Bonnet" (plus the ⭙ glyph)
 *  is the second channel, so the signal survives grayscale and screen readers.
 *
 *  Voice (E): uppercase bold sans label is chrome/label voice, correct.
 *
 *  Sizing (F): the mockup's ~9px is below the text-xs (~10.5px) floor; we use
 *  text-xs (+1.5px, on the fixed scale) rather than mint a sub-xs token.
 *  The dynamic color-mix lives in an inline `style` using token CSS vars
 *  (mirroring PhaseChip's dynamic phase color), never a raw hex/oklch literal. */
export function BonnetBadge({ className = "" }: BonnetBadgeProps): JSX.Element {
  return (
    <span
      data-testid="bonnet-badge"
      className={cx(
        "inline-flex items-center gap-0.5 text-xs font-bold uppercase tracking-wide rounded-full px-1.5 py-px",
        className,
      )}
      style={{
        color: "var(--color-accent)",
        borderWidth: 1,
        borderStyle: "solid",
        borderColor: "color-mix(in oklab, var(--color-accent) 40%, var(--color-hair))",
        background: "color-mix(in oklab, var(--color-accent) 8%, transparent)",
      }}
    >
      <span aria-hidden="true">⭙</span>
      Bonnet
    </span>
  );
}
