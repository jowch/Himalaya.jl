interface ToggleSwitchProps {
  checked: boolean;
  onChange: (v: boolean) => void;
  label: string;
  /** Visually hide the label text while keeping it as the accessible name +
   *  screen-reader text, for dense standalone use (e.g. a toolbar). */
  hideLabel?: boolean;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A 32x18 pill switch with a 14px plate knob that slides 14px when on, plus a
 *  text label beside it. Built pseudo-free (a `<button role="switch">` + spans)
 *  so it's testable and fully tokenized.
 *
 *  A — the on/off state carries via knob POSITION (`translate-x-3.5`) plus
 *  `role="switch"`/`aria-checked`, never hue alone; it reads in grayscale. C —
 *  focus is the focus-visible accent outline on the switch button (hover N/A on
 *  the track). D — the 32x18 control is small, but the FULL `<label>` row
 *  (switch + text) is the hit target via the wrapping label + `cursor-pointer`.
 *  G — only `transition-colors` (track) + `transition-transform` (knob), both
 *  allowed; no layout animation. H — the knob is FLAT: it uses a
 *  `border border-hair-strong` hairline ring for definition, NOT a `shadow-*`
 *  drop shadow. Rationale: checklist H reserves shadow for the Card plate only;
 *  a hairline ring keeps the knob flat AND matches the Slider thumb's 1px
 *  hair-strong ring. */
export function ToggleSwitch({
  checked,
  onChange,
  label,
  hideLabel = false,
  className,
}: ToggleSwitchProps): JSX.Element {
  return (
    <label
      data-testid="toggle-switch"
      className={cx("inline-flex items-center gap-2.5 cursor-pointer", className)}
    >
      <button
        type="button"
        role="switch"
        aria-checked={checked}
        aria-label={label}
        onClick={() => onChange(!checked)}
        className={cx(
          "relative w-8 h-[18px] rounded-full transition-colors flex-shrink-0",
          "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
          checked ? "bg-accent" : "bg-hair-strong",
        )}
      >
        <span
          className={cx(
            "absolute top-0.5 left-0.5 w-3.5 h-3.5 rounded-full bg-plate border border-hair-strong transition-transform",
            checked && "translate-x-3.5",
          )}
        />
      </button>
      <span
        data-role="switch-label"
        data-hidden={hideLabel ? "true" : undefined}
        className={cx("text-base font-semibold text-ink", hideLabel && "sr-only")}
      >
        {label}
      </span>
    </label>
  );
}
