import { Slider, Input } from "../ui";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface LatticeParamControlProps {
  value: string;
  min: number;
  max: number;
  step?: number;
  onValueChange: (v: string) => void;
  unit?: string;
  /** a11y name for the slider + number field. */
  label?: string;
  /** Second-channel an out-of-range / non-finite value: flags the number field
   *  `aria-invalid` (red border) so a disabled Add button is explained at the
   *  field, mirroring the trace q-add field. */
  invalid?: boolean;
  className?: string;
}

/** A bare range slider (fills) + a mono number field (64px) + a unit label.
 *  Presentational: `value` is a string, `onValueChange` returns a string.
 *  The slider emits a number which is stringified; the Input emits a string
 *  directly. No visible label row — a11y name is carried via `ariaLabel`. */
export function LatticeParamControl({
  value,
  min,
  max,
  step,
  onValueChange,
  unit = "Å",
  label = "lattice parameter",
  invalid = false,
  className,
}: LatticeParamControlProps): JSX.Element {
  return (
    <div data-testid="lattice-param" className={cx("flex items-center gap-2.5 flex-1 min-w-0", className)}>
      <Slider
        value={Number(value)}
        min={min}
        max={max}
        {...(step !== undefined ? { step } : {})}
        onChange={(n) => onValueChange(String(n))}
        ariaLabel={label}
        className="flex-1 min-w-0"
      />
      <Input
        testId="lattice-num"
        type="number"
        mono
        inputSize="sm"
        value={value}
        onValueChange={onValueChange}
        min={min}
        max={max}
        invalid={invalid}
        {...(step !== undefined ? { step } : {})}
        aria-label={`${label} value`}
        className="w-16 shrink-0"
      />
      <span className="text-caption text-ink-soft font-mono shrink-0">{unit}</span>
    </div>
  );
}
