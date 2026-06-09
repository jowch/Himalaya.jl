import { useRef } from "react";
import type { ReactNode, KeyboardEvent } from "react";

/**
 * SegmentedControl<T> — the canonical single-select button group.
 *
 * Closed look / open placement (the Button.tsx pattern): the library owns
 * appearance via `variant`/`size` and the canonical ink-on-paper active fill
 * (DESIGN.md §211/§240 — the active segment is `bg-ink text-paper`, never the
 * recessed `bg-paper-sunk` "L-5/B-A defect"). The consumer's `className` is
 * PLACEMENT-ONLY (margin / width / flex position); the design guard bans
 * appearance utilities there.
 *
 * `role` drives child semantics: "group" -> button[aria-pressed];
 * "radiogroup" -> button[role=radio][aria-checked] with WAI-ARIA roving
 * tabindex + arrow-key navigation. The container always reflects the active
 * value via `data-state` (and, if `stateAttr` is given, an extra aliased
 * attribute e.g. `data-mode` for GroupingModeToggle's E2E contract); each
 * segment reflects `data-active`/`data-value`.
 *
 * 44px touch target is folded in via the `.seg-segment` class behind a
 * `@media (pointer:coarse)` rule in styles.css, so dense desktop toolbars
 * keep their compact density on a fine pointer.
 */

export interface SegmentOption<T extends string> {
  value: T;
  label: ReactNode;
  /** Disabled individual segment (e.g. Loupe with no sample selected). */
  disabled?: boolean;
  /** Per-segment title/tooltip. */
  title?: string;
  /** Stable test id for E2E, applied to the segment button. */
  testId?: string;
}

export type SegmentedVariant = "bordered" | "plain";
export type SegmentedSize = "xs" | "sm" | "md";

export interface SegmentedControlProps<T extends string> {
  options: ReadonlyArray<SegmentOption<T>>;
  value: T;
  onChange: (next: T) => void;
  /** group | radiogroup — drives child role + keyboard model. Default "group". */
  role?: "group" | "radiogroup";
  variant?: SegmentedVariant;
  size?: SegmentedSize;
  /** Required: names the control for assistive tech. */
  "aria-label": string;
  /** Container test id. */
  testId?: string;
  /**
   * Extra container attribute name aliasing the active value, in addition to
   * the always-present `data-state` (e.g. "data-mode" for GroupingModeToggle).
   */
  stateAttr?: string;
  /** Makes the container full-width and each segment flex-1 (equal-width fill). */
  stretch?: boolean;
  /** PLACEMENT-ONLY: margin / width / flex-grid position. No appearance utils. */
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

const containerClass: Record<SegmentedVariant, string> = {
  bordered: "inline-flex overflow-hidden rounded border border-hair-strong divide-x divide-hair-strong",
  plain: "inline-flex items-center gap-1",
};

// `xs` is a known-dense in-panel toggle (focus-plot mini comb switch). Its
// effective hit area is enlarged by the toolbar row it sits in; the segment
// itself is intentionally below the 44px standalone touch target. `sm`/`md`
// remain the default comfortable sizes.
const sizeClass: Record<SegmentedSize, string> = {
  xs: "text-xs px-2 py-1",
  sm: "text-xs px-3 py-1.5",
  md: "text-sm px-3.5 py-2",
};

const segmentBaseClass: Record<SegmentedVariant, string> = {
  bordered: "",
  plain: "rounded border border-transparent",
};

export function SegmentedControl<T extends string>({
  options,
  value,
  onChange,
  role = "group",
  variant = "bordered",
  size = "sm",
  "aria-label": ariaLabel,
  testId,
  stateAttr,
  stretch = false,
  className = "",
}: SegmentedControlProps<T>): JSX.Element {
  const isRadio = role === "radiogroup";
  // Refs to the segment buttons so radiogroup arrow-nav can move DOM focus to
  // the newly-selected radio (WAI-ARIA: focus follows selection). Without this,
  // activeElement is left on a now-tabindex=-1 segment after an arrow key.
  const segRefs = useRef<Array<HTMLButtonElement | null>>([]);

  const move = (delta: number, from: number): void => {
    const n = options.length;
    if (n === 0) return;
    let i = from;
    for (let step = 0; step < n; step++) {
      i = (i + delta + n) % n;
      if (!options[i].disabled) {
        onChange(options[i].value);
        segRefs.current[i]?.focus();
        return;
      }
    }
  };

  const onKeyDown = (e: KeyboardEvent<HTMLButtonElement>, idx: number): void => {
    if (!isRadio) return;
    if (e.key === "ArrowRight" || e.key === "ArrowDown") {
      e.preventDefault();
      move(1, idx);
    } else if (e.key === "ArrowLeft" || e.key === "ArrowUp") {
      e.preventDefault();
      move(-1, idx);
    }
  };

  const containerProps: Record<string, string> = {
    "data-state": value,
    "data-variant": variant,
    "data-size": size,
  };
  if (stateAttr) containerProps[stateAttr] = value;
  if (testId) containerProps["data-testid"] = testId;

  return (
    <div
      role={role}
      aria-label={ariaLabel}
      className={cx(containerClass[variant], stretch && "flex w-full", className)}
      {...containerProps}
    >
      {options.map((opt, idx) => {
        const active = opt.value === value;
        const selectedProps = isRadio
          ? { role: "radio" as const, "aria-checked": active, tabIndex: active ? 0 : -1 }
          : { "aria-pressed": active };
        return (
          <button
            key={opt.value}
            ref={(el) => {
              segRefs.current[idx] = el;
            }}
            type="button"
            disabled={opt.disabled}
            title={opt.title}
            data-testid={opt.testId}
            data-value={opt.value}
            data-active={active ? "true" : "false"}
            onClick={() => onChange(opt.value)}
            onKeyDown={(e) => onKeyDown(e, idx)}
            className={cx(
              "seg-segment transition-colors focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
              sizeClass[size],
              segmentBaseClass[variant],
              active
                ? "bg-ink text-paper"
                : "text-ink-soft hover:text-ink hover:bg-paper-sunk",
              opt.disabled && "opacity-50 cursor-not-allowed",
              stretch && "flex-1",
            )}
            {...selectedProps}
          >
            {opt.label}
          </button>
        );
      })}
    </div>
  );
}
