import { useEffect, useRef } from "react";
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
 * Semantics: `role="menu"` container. Each option is a `role="menuitem"`
 * (action menus: export, nav) OR a `role="menuitemradio"` with `aria-checked`
 * (VALUE-SELECTOR menus — those given an `activeValue`, e.g. the ordering
 * dropdown). The radio role tells assistive tech WHICH option is the current
 * value; a plain menuitem cannot. The `role="menu"` container validly holds
 * either. Keyboard model mirrors SegmentedControl's radiogroup roving focus —
 * refs array + `.focus()` on ArrowUp/Down — plus Escape -> onClose. A click
 * selects then closes. Renders nothing when `!open`.
 *
 * APG menu-button: opening the menu moves DOM focus INTO it — to the active
 * item when `activeValue` names an enabled option, else the first enabled
 * item; `initialFocus="last"` targets the last enabled item instead, so owners
 * can implement ArrowUp-on-trigger-opens-focusing-last. Mouse clicks on the
 * trigger land focus in the menu too (the APG makes no pointer/keyboard
 * distinction on open); `focus-visible` keeps the focus ring keyboard-only, so
 * mouse users see no ring flash. Returning focus to the TRIGGER on close is
 * the owner's job (Menu has no trigger ref) — see Field / ExportButton.
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
  /** The currently-active option's value (gets `data-active` + initial focus). */
  activeValue?: T;
  /** Where focus lands on open (APG menu-button). Default/"first": the active
   *  item if enabled, else the first enabled item. "last": the last enabled
   *  item — owners pass this when ArrowUp opened the menu. */
  initialFocus?: "first" | "last";
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
  initialFocus,
  className = "",
}: MenuProps<T>): JSX.Element | null {
  // Refs to the menuitem buttons so arrow-nav can move DOM focus among them,
  // mirroring SegmentedControl's radiogroup roving-focus idiom.
  const itemRefs = useRef<Array<HTMLButtonElement | null>>([]);

  // APG menu-button: when the menu opens, move DOM focus to a menuitem. The
  // menu mounts on open (we render null when closed), so an effect keyed on
  // `open` fires exactly on the open transition, after the item refs exist.
  // Deliberately NOT keyed on options/activeValue/initialFocus — focus moves
  // once per open, never on later re-renders of an already-open menu.
  useEffect(() => {
    if (!open) return;
    let target = -1;
    if (initialFocus === "last") {
      for (let i = options.length - 1; i >= 0; i--) {
        if (!options[i].disabled) {
          target = i;
          break;
        }
      }
    } else {
      // Active item wins (so the menu opens "on" the current value); a
      // disabled active option falls through to first-enabled.
      target = options.findIndex((o) => o.value === activeValue && !o.disabled);
      if (target < 0) target = options.findIndex((o) => !o.disabled);
    }
    if (target >= 0) itemRefs.current[target]?.focus();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [open]);

  if (!open) return null;

  // A value-selector menu (the consumer passes an activeValue) marks the
  // current option with role=menuitemradio + aria-checked; an action menu
  // (no activeValue: export, nav) keeps the bare menuitem with no checked state.
  const isValueSelector = activeValue !== undefined;

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
      // APG: Escape dismisses only the innermost popup. Stop the event here
      // so a host modal's document-level Escape listener (ModalShell) does
      // not double-close.
      e.stopPropagation();
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
            role={isValueSelector ? "menuitemradio" : "menuitem"}
            type="button"
            // UI-MENUTAB — the menu is a single tab stop (APG menu pattern):
            // items are NOT in the tab order, so one Tab from a focused item
            // leaves the whole popup instead of walking item-to-item, which
            // lets the owner's focus-out handler close it (APG: Tab closes the
            // menu). Focus still moves here via the open effect + arrow nav
            // (imperative .focus(), unaffected by tabIndex).
            tabIndex={-1}
            disabled={o.disabled}
            {...(isValueSelector ? { "aria-checked": active } : {})}
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
