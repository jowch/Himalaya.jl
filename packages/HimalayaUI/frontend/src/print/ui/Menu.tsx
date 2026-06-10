import { useRef } from "react";
import type { ReactNode, KeyboardEvent } from "react";

/**
 * Menu<T> — the plate dropdown popover (closed look / open placement).
 *
 * Renders as the `.card` Plate-Lift surface (DESIGN.md §H — Menu is plate-like,
 * so the lifted plate is the correct read here): an absolutely-positioned
 * popover that the consumer anchors via a `relative` parent. The library owns
 * appearance (the `.card` plate + menuitem hover/focus/active recipe); the
 * consumer's `className` is PLACEMENT-ONLY and lands LAST.
 *
 * Semantics: `role="menu"` with `role="menuitem"` children. Keyboard model
 * mirrors SegmentedControl's radiogroup roving focus — refs array + `.focus()`
 * on ArrowUp/Down — plus Escape -> onClose. A click selects then closes.
 * Renders nothing when `!open`.
 */

export interface MenuOption<T extends string> {
  value: T;
  label: ReactNode;
  /** Disabled individual option (not selectable, skipped by arrow nav). */
  disabled?: boolean;
}

export interface MenuProps<T extends string> {
  open: boolean;
  options: ReadonlyArray<MenuOption<T>>;
  onSelect: (value: T) => void;
  onClose: () => void;
  /** Required: names the menu for assistive tech. */
  "aria-label": string;
  /** The currently-active option's value (gets `data-active`). */
  activeValue?: T;
  /** PLACEMENT-ONLY: margin / width / position. No appearance utils. */
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export function Menu<T extends string>({
  open,
  options,
  onSelect,
  onClose,
  "aria-label": ariaLabel,
  activeValue,
  className = "",
}: MenuProps<T>): JSX.Element | null {
  // Refs to the menuitem buttons so arrow-nav can move DOM focus among them,
  // mirroring SegmentedControl's radiogroup roving-focus idiom.
  const itemRefs = useRef<Array<HTMLButtonElement | null>>([]);

  if (!open) return null;

  const move = (delta: number, from: number): void => {
    const n = options.length;
    if (n === 0) return;
    let i = from;
    for (let step = 0; step < n; step++) {
      i = (i + delta + n) % n;
      if (!options[i].disabled) {
        itemRefs.current[i]?.focus();
        return;
      }
    }
  };

  const onItemKeyDown = (e: KeyboardEvent<HTMLButtonElement>, idx: number): void => {
    if (e.key === "ArrowDown") {
      e.preventDefault();
      move(1, idx);
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      move(-1, idx);
    }
  };

  const onMenuKeyDown = (e: KeyboardEvent<HTMLDivElement>): void => {
    if (e.key === "Escape") {
      e.preventDefault();
      onClose();
    }
  };

  return (
    <div
      role="menu"
      aria-label={ariaLabel}
      data-testid="menu"
      onKeyDown={onMenuKeyDown}
      className={cx("card absolute z-20 mt-1 min-w-[180px] py-1", className)}
    >
      {options.map((o, idx) => {
        const active = o.value === activeValue;
        return (
          <button
            key={o.value}
            ref={(el) => {
              itemRefs.current[idx] = el;
            }}
            role="menuitem"
            type="button"
            disabled={o.disabled}
            data-value={o.value}
            data-active={active ? "true" : undefined}
            onClick={() => {
              onSelect(o.value);
              onClose();
            }}
            onKeyDown={(e) => onItemKeyDown(e, idx)}
            className={cx(
              "flex w-full items-center px-3 py-1.5 text-sm text-left transition-colors",
              active
                ? "text-ink bg-paper-sunk"
                : "text-ink-soft hover:text-ink hover:bg-paper-sunk",
              "disabled:opacity-40 disabled:cursor-not-allowed focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
            )}
          >
            {o.label}
          </button>
        );
      })}
    </div>
  );
}
