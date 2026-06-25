import { cx } from "../../lib/cx";
interface SignalBarsProps {
  /** Strength reading, on the same scale as `max`. Mapped onto 5 bars. */
  value: number;
  /** Top of the scale `value` is read against. Default 5. */
  max?: number;
  /** PLACEMENT ONLY — margin / inline position. Appearance utilities are
   *  banned by the lint:design guard; the look lives in this file. */
  className?: string;
}


/** A 5-bar signal-strength indicator. The active COUNT = round(value/max * 5),
 *  clamped to 0..5.
 *
 *  Second-channel (checklist A): strength is carried by the COUNT/position of
 *  filled bars (more bars = stronger), never by hue. The fill is a tonal step —
 *  `ink-soft` (on) vs `hair-strong` (off) — so it reads in grayscale; the
 *  `aria-label` ("Signal N of 5") restates the count for screen readers. There
 *  is no chromatic code.
 *
 *  Interaction (C/D): N/A — a non-interactive read-only indicator. No hover/
 *  focus/active/disabled states and no touch target (it carries no action).
 *
 *  Flat (H): tonal bars + no shadow; not a plate.
 *
 *  Sizing (F): `w-[5px] h-[11px] rounded-[1px]` are fixed semantic bar
 *  dimensions (the mockup's `.signal-bars` geometry), allowed in print/ui. */
export function SignalBars({ value, max, className = "" }: SignalBarsProps): JSX.Element {
  const denom = max ?? 5;
  const on = Math.max(0, Math.min(5, Math.round((value / denom) * 5)));
  return (
    <span
      data-testid="signal-bars"
      role="img"
      aria-label={`Signal ${on} of 5`}
      className={cx("inline-flex gap-0.5 items-end", className)}
    >
      {Array.from({ length: 5 }, (_, i) => (
        <i
          key={i}
          data-on={i < on ? "true" : "false"}
          // Fixed semantic bar dimensions (mockup .signal-bars geometry).
          className={cx("w-[5px] h-[11px] rounded-[1px]", i < on ? "bg-ink-soft" : "bg-hair-strong")}
        />
      ))}
    </span>
  );
}
