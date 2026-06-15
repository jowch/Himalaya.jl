import type { ReactNode } from "react";

interface SliderProps {
  value: number;
  min: number;
  max: number;
  step?: number;
  onChange: (v: number) => void;
  label?: string;
  valueDisplay?: ReactNode;
  /** Sets `aria-label` on the `<input>` without rendering a visible label row.
   *  Useful for bare sliders embedded in a composite that provides its own
   *  visible label. Precedence: `ariaLabel` wins over `label` for the input's
   *  `aria-label`; the visible row is still gated on `label || valueDisplay`. */
  ariaLabel?: string;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A styled `<input type="range">` with a thin hair-strong track and an accent
 *  thumb, plus an optional label/value row above it. The custom thumb + track
 *  require real CSS pseudo-elements (`::-webkit-slider-thumb`,
 *  `::-moz-range-thumb`) that Tailwind cannot express, so appearance lives in
 *  the `.print-range` class in `styles.css` (the legitimate home for
 *  pseudo-element styling). Consumer `className` is placement-only and lands
 *  last on the outer column.
 *
 *  A — the slider's value carries via the thumb POSITION plus the optional mono
 *  `valueDisplay`, never hue. B — the accent thumb is the live control (the
 *  rationed, semantic accent). C — focus is the focus-visible thumb outline
 *  (authored in the CSS); hover N/A on a range thumb. G — no motion. */
export function Slider({
  value,
  min,
  max,
  step,
  onChange,
  label,
  valueDisplay,
  ariaLabel,
  className,
}: SliderProps): JSX.Element {
  const a11yName = ariaLabel ?? label;
  return (
    <div className={cx("flex flex-col gap-1", className)}>
      {(label || valueDisplay) && (
        <div className="flex justify-between items-baseline text-xs text-ink-soft">
          {label && <span>{label}</span>}
          {valueDisplay && <span className="font-mono text-ink-soft">{valueDisplay}</span>}
        </div>
      )}
      <input
        type="range"
        data-testid="slider"
        className="print-range"
        value={value}
        min={min}
        max={max}
        {...(step !== undefined ? { step } : {})}
        onChange={(e) => onChange(Number(e.target.value))}
        {...(a11yName ? { "aria-label": a11yName } : {})}
      />
    </div>
  );
}
