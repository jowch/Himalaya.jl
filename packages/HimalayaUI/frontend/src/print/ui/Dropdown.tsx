import { useRef, useState } from "react";
import type { ReactNode } from "react";
import { Menu } from "./Menu";
import type { MenuOption } from "./Menu";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface DropdownProps<T extends string> {
  /** Required: names the control for assistive tech. */
  "aria-label": string;
  options: ReadonlyArray<MenuOption<T>>;
  /** The current value. `undefined` shows the placeholder. */
  value: T | undefined;
  onChange: (value: T) => void;
  /** Shown on the trigger when no option matches `value`. */
  placeholder?: ReactNode;
  /** Disables the trigger. */
  disabled?: boolean;
  /** PLACEMENT-ONLY: margin / width / position. No appearance utils. */
  className?: string;
}

/**
 * Dropdown<T> — a labelled select-trigger that opens a {@link Menu} value-
 * selector popover (closed look / open placement). The trigger shows the active
 * option's label (or `placeholder`) + a chevron; clicking toggles the Menu,
 * which owns the list look + roving focus + Escape. Selecting closes the menu
 * and returns focus to the trigger (Menu leaves focus-return to the owner).
 *
 * The reusable replacement for ad-hoc `<select>` + `Menu` wirings (spec: "turn
 * the dropdown into a component and reuse it"). Consumed by DirectoryPickerField
 * here and the grouping/sources menus in E2.
 */
export function Dropdown<T extends string>({
  "aria-label": ariaLabel,
  options,
  value,
  onChange,
  placeholder,
  disabled = false,
  className = "",
}: DropdownProps<T>): JSX.Element {
  const [open, setOpen] = useState(false);
  const triggerRef = useRef<HTMLButtonElement>(null);

  const active = options.find((o) => o.value === value);
  const triggerLabel = active ? active.label : (placeholder ?? "");

  const close = (): void => {
    setOpen(false);
    triggerRef.current?.focus();
  };

  return (
    <div className={cx("relative inline-flex", className)}>
      <button
        ref={triggerRef}
        type="button"
        data-testid="dropdown-trigger"
        aria-haspopup="menu"
        aria-expanded={open}
        aria-label={ariaLabel}
        disabled={disabled}
        onClick={() => setOpen((v) => !v)}
        className={cx(
          "inline-flex items-center justify-between gap-2 min-w-0",
          "rounded-sm border border-hair-strong bg-plate px-2.5 py-1.5",
          "text-sm text-ink transition-colors hover:border-accent",
          "focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2",
          "disabled:opacity-40 disabled:cursor-not-allowed",
        )}
      >
        <span className={cx("truncate", active ? "text-ink" : "text-ink-soft")}>
          {triggerLabel}
        </span>
        <span aria-hidden="true" className="text-ink-faint">▾</span>
      </button>
      <Menu<T>
        open={open}
        options={options}
        {...(value !== undefined ? { activeValue: value } : {})}
        onSelect={onChange}
        onClose={close}
        aria-label={ariaLabel}
        className="left-0 top-full w-full"
      />
    </div>
  );
}
