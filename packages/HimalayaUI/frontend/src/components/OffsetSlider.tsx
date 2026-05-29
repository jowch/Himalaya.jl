interface OffsetSliderProps {
  /** Current offset multiplier (0.4..1.4). */
  value: number;
  onChange: (next: number) => void;
  /** Greyed + inert when the offset has no effect (e.g. heatmap mode). */
  disabled?: boolean;
}

/**
 * OffsetSlider — the series-builder "trace offset" control (R8 / B-F). Scales
 * the vertical spacing of the waterfall stack. The page maps this 0.4..1.4
 * range to a working-band fraction via `offsetToBandFraction`. Mockup:
 * `series-builder.html` `.slider-row` / `#offset`.
 */
export function OffsetSlider({ value, onChange, disabled = false }: OffsetSliderProps): JSX.Element {
  return (
    <div
      className={`flex flex-col gap-1.5 ${disabled ? "opacity-40" : ""}`}
      data-testid="offset-slider-row"
    >
      <div className="flex items-baseline justify-between">
        <span className="text-xs font-semibold text-ink">Trace offset</span>
        <span className="text-xs font-semibold text-ink" data-testid="offset-slider-value">
          {value.toFixed(2)}&times;
        </span>
      </div>
      <input
        type="range"
        data-testid="offset-slider"
        aria-label="Trace offset"
        min="0.4"
        max="1.4"
        step="0.05"
        value={value}
        disabled={disabled}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="h-[3px] w-full appearance-none rounded-full bg-hair-strong accent-print-accent disabled:cursor-not-allowed cursor-pointer"
      />
    </div>
  );
}
