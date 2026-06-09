import { useEffect, useRef, useState } from "react";
import { Menu } from "./Menu";
import type { MenuOption } from "./Menu";

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
  /** data-testid override (default "field"). */
  testId?: string;
  className?: string; // PLACEMENT ONLY
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** The `.field` ordering-variable control: a bordered, clickable row showing a
 *  value on the left and a `▾` chevron (ink-faint) on the right.
 *
 *  In DROPDOWN mode (`options` given) the field is self-contained: it owns its
 *  open state, renders the `Menu` primitive anchored to a `relative` wrapper,
 *  and closes on select / outside-click / Escape (Menu owns Escape + arrow nav).
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
  testId,
  className,
}: FieldProps): JSX.Element {
  const showPlaceholder = value === "" && placeholder != null;
  const isDropdown = options != null;
  const [open, setOpen] = useState(false);
  const wrapRef = useRef<HTMLSpanElement | null>(null);

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
        <span className={showPlaceholder ? "text-ink-soft" : "text-ink"}>
          {showPlaceholder ? placeholder : value}
        </span>
      </div>
    );
  }

  const trigger = (
    <button
      type="button"
      data-testid={testId ?? "field"}
      {...(isDropdown
        ? { "aria-haspopup": "menu" as const, "aria-expanded": open }
        : {})}
      onClick={() => {
        if (isDropdown) setOpen((o) => !o);
        onClick?.();
      }}
      className={cx(
        "w-full flex items-center justify-between border border-hair-strong bg-plate rounded px-3 py-2 text-meta font-semibold cursor-pointer",
        isDropdown ? undefined : className,
      )}
    >
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
    <span ref={wrapRef} className={cx("relative inline-block w-full", className)}>
      {trigger}
      <Menu<string>
        open={open}
        options={options}
        onSelect={(v) => {
          onSelect?.(v);
          setOpen(false);
        }}
        onClose={() => setOpen(false)}
        aria-label={menuLabel ?? "Choose an option"}
        activeValue={value}
        className="w-full"
      />
    </span>
  );
}
