import { useEffect, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Menu } from "./Menu";
import type { MenuOption } from "./Menu";
import { cx } from "../../lib/cx";

export interface FieldProps {
  value: string;
  /** When provided → self-contained dropdown (internal open + Menu). */
  options?: ReadonlyArray<MenuOption<string>>;
  /** Selection handler (dropdown mode). */
  onSelect?: (value: string) => void;
  /** Trigger handler (no-options fallback mode). */
  onClick?: () => void;
  placeholder?: string;
  /** a11y name for the menu (dropdown mode). Default "Choose an option". */
  menuLabel?: string;
  /** Visually-hidden descriptive label rendered before the value, so the
   *  control's accessible name reads "{srLabel} {value}" (WCAG 4.1.2).
   *  Without it the name is just the bare value. */
  srLabel?: string;
  /** data-testid override (default "field"). */
  testId?: string;
  className?: string; // PLACEMENT ONLY
}


/** The `.field` ordering-variable control: a bordered, clickable row showing a
 *  value on the left and a `▾` chevron (ink-faint) on the right.
 *
 *  In DROPDOWN mode (`options` given) the field is self-contained: it owns its
 *  open state, renders the `Menu` primitive anchored to a `relative` wrapper,
 *  and closes on select / outside-click / Escape (Menu owns in-menu Escape +
 *  arrow nav; Menu's mount effect moves focus into the menu on open).
 *
 *  APG menu-button keyboard pattern (the owner half):
 *  - ArrowDown on the closed trigger opens + focuses the active/first item;
 *    ArrowUp opens + focuses the LAST item (both preventDefault — no scroll).
 *  - Escape on the trigger while the menu is open closes it (focus is already
 *    on the trigger, so nothing moves).
 *  - Any keyboard-reachable close (Escape inside the menu, selecting an item)
 *    RETURNS focus to the trigger. Outside-pointerdown close deliberately does
 *    NOT — the pointer user clicked elsewhere on purpose; yanking focus back
 *    to the trigger would fight their intent.
 *
 *  In TRIGGER mode (no `options`) it is the bare button whose click fires the
 *  passed `onClick` — the menu is owned elsewhere (today's behaviour).
 *
 *  When `value` is empty and a `placeholder` is given, the placeholder reads in
 *  ink-faint. */
export function Field({
  value,
  options,
  onSelect,
  onClick,
  placeholder,
  menuLabel,
  srLabel,
  testId,
  className,
}: FieldProps): JSX.Element {
  const showPlaceholder = value === "" && placeholder != null;
  const isDropdown = options != null;
  const [open, setOpen] = useState(false);
  // Where Menu's mount effect puts focus: "first" = active-else-first (click,
  // ArrowDown); "last" = last enabled item (ArrowUp). APG menu-button.
  const [menuFocus, setMenuFocus] = useState<"first" | "last">("first");
  const wrapRef = useRef<HTMLSpanElement | null>(null);
  const triggerRef = useRef<HTMLButtonElement | null>(null);

  // Outside-click closes the menu (mirrors Popover's outside-pointerdown
  // pattern); bound only while open.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: PointerEvent): void => {
      if (wrapRef.current && !wrapRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, [open]);

  // controls-don't-lie: a Field with neither `options` (dropdown) nor `onClick`
  // (trigger) is INERT — it must read as a static read-only value, not as a
  // clickable dropdown. Render a plain non-interactive row: no `▾` chevron, no
  // pointer cursor, not a `<button>`. The chevron + pointer are reserved for the
  // genuinely-actionable dropdown / trigger forms below.
  const isStatic = !isDropdown && onClick === undefined;
  if (isStatic) {
    return (
      <div
        data-testid={testId ?? "field"}
        className={cx(
          "w-full flex items-center justify-between border border-hair-strong bg-plate rounded px-3 py-2 text-meta font-semibold",
          className,
        )}
      >
        {srLabel != null && <span className="sr-only">{srLabel} </span>}
        <span className={showPlaceholder ? "text-ink-soft" : "text-ink"}>
          {showPlaceholder ? placeholder : value}
        </span>
      </div>
    );
  }

  // Close + return focus to the trigger — every keyboard-reachable close path
  // (Escape inside the menu, item select) routes through here so the user is
  // never dropped at <body>. Outside-pointerdown bypasses it on purpose.
  const closeAndRefocus = (): void => {
    setOpen(false);
    triggerRef.current?.focus();
  };

  // APG menu-button trigger keys (dropdown mode only):
  // ArrowDown opens focusing active/first; ArrowUp opens focusing last;
  // Escape closes an open menu whose focus is still on the trigger (without
  // this, Escape here was inert — Menu only hears keys once focus is inside).
  // Enter/Space already open via the native button click.
  const onTriggerKeyDown = (e: KeyboardEvent<HTMLButtonElement>): void => {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault(); // arrows must not scroll the page
      setMenuFocus(e.key === "ArrowDown" ? "first" : "last");
      setOpen(true);
    } else if (e.key === "Escape" && open) {
      e.preventDefault();
      e.stopPropagation(); // innermost popup only — don't double-close a host modal
      setOpen(false); // focus is already on the trigger — nothing to restore
    }
  };

  const trigger = (
    <button
      ref={triggerRef}
      type="button"
      data-testid={testId ?? "field"}
      {...(isDropdown
        ? {
            "aria-haspopup": "menu" as const,
            "aria-expanded": open,
            onKeyDown: onTriggerKeyDown,
          }
        : {})}
      onClick={() => {
        if (isDropdown) {
          setMenuFocus("first"); // click-open focuses active-else-first
          setOpen((o) => !o);
        }
        onClick?.();
      }}
      className={cx(
        "w-full flex items-center justify-between border border-hair-strong bg-plate rounded px-3 py-2 text-meta font-semibold cursor-pointer",
        isDropdown ? undefined : className,
      )}
    >
      {srLabel != null && <span className="sr-only">{srLabel} </span>}
      <span className={showPlaceholder ? "text-ink-soft" : "text-ink"}>
        {showPlaceholder ? placeholder : value}
      </span>
      <span className="text-ink-faint" aria-hidden="true">
        ▾
      </span>
    </button>
  );

  if (!isDropdown) return trigger;

  return (
    <span
      ref={wrapRef}
      // UI-MENUTAB — APG menu-button: Tab (Shift+Tab) closes the menu. Mirror
      // the outside-pointerdown close on the FOCUS axis — when focus leaves the
      // whole widget (relatedTarget outside wrapRef, or null) close it. Focus
      // moving to the trigger or between items stays inside wrapRef, so a
      // toggle-click never closes-then-reopens. Plain setOpen (no refocus): the
      // user is already on their way out.
      onBlur={(e) => {
        if (wrapRef.current && !wrapRef.current.contains(e.relatedTarget as Node | null)) {
          setOpen(false);
        }
      }}
      className={cx("relative inline-block w-full", className)}
    >
      {trigger}
      <Menu<string>
        open={open}
        options={options}
        onSelect={(v) => onSelect?.(v)}
        onClose={closeAndRefocus}
        aria-label={menuLabel ?? "Choose an option"}
        activeValue={value}
        initialFocus={menuFocus}
        className="w-full"
      />
    </span>
  );
}
