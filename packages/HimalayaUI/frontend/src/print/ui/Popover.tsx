import { cloneElement, useEffect, useId, useRef, useState } from "react";
import type { KeyboardEvent, ReactElement, ReactNode, RefObject } from "react";
import { cx } from "../../lib/cx";

export interface PopoverProps {
  /** The trigger element. Cloned to wire click + aria. Must be a single
   *  focusable element (a `<button>`). */
  trigger: ReactElement;
  /** Popover body content (may contain interactive elements). */
  children: ReactNode;
  /** Side preference. Default `"bottom"`. */
  side?: "top" | "bottom";
  /** Accessible name for the popover dialog. */
  label?: string;
  /** Element to focus when the popover opens, instead of the panel itself.
   *  For a search-first popover (e.g. a combobox), point this at the inner
   *  input so keyboard users can type immediately. The element mounts with the
   *  panel, so its ref is set by the time the open effect runs. Falls back to
   *  the panel when absent (the default). */
  initialFocusRef?: RefObject<HTMLElement | null>;
  /** Stretch the trigger wrapper to the full width of its container (block flex
   *  instead of the default inline-flex). For a full-width combobox bar; the
   *  trigger itself must also be `w-full`. Default false (inline). */
  fullWidth?: boolean;
  /** PLACEMENT-ONLY on the popover panel. */
  className?: string;
}


/**
 * Popover — a click/focus-triggered LIGHT-plate popover holding ARBITRARY,
 * possibly interactive content. The W3C ARIA counterpart to {@link Tooltip}:
 * where a hover tooltip may hold only non-interactive text, this dialog can
 * carry focusable controls (the whole point), so it opens on click and stays
 * open until Escape / outside-pointerdown / re-click.
 *
 * Surface mirrors {@link Menu} — the `.card` Plate-Lift plate, absolutely
 * positioned, anchored by a `relative` wrapper, with the `.anim-pal-scale`
 * popover entry. It is NOT the dark Tooltip exception.
 *
 * A11y: the trigger is cloned to carry `aria-haspopup="dialog"`,
 * `aria-expanded`, and `aria-controls`; the panel is `role="dialog"` with
 * `aria-label`. Escape closes and returns focus to the trigger. Ids come from
 * React `useId()` — no `Math.random()`/`Date.now()`.
 */
export function Popover({
  trigger,
  children,
  side = "bottom",
  label,
  initialFocusRef,
  fullWidth = false,
  className = "",
}: PopoverProps): JSX.Element {
  const panelId = useId();
  const [open, setOpen] = useState(false);
  const wrapperRef = useRef<HTMLSpanElement>(null);
  const triggerRef = useRef<HTMLElement | null>(null);
  const panelRef = useRef<HTMLDivElement>(null);

  // Close on outside pointerdown. Bound only while open.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: PointerEvent): void => {
      if (!wrapperRef.current) return;
      if (!wrapperRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, [open]);

  // Move focus into the panel when it opens, so keyboard users land on the
  // revealed content (and Escape has a focused target). A search-first popover
  // can redirect that focus to its input via `initialFocusRef`.
  useEffect(() => {
    if (open) (initialFocusRef?.current ?? panelRef.current)?.focus();
  }, [open, initialFocusRef]);

  const close = (): void => {
    setOpen(false);
    triggerRef.current?.focus();
  };

  const existingRef = (
    trigger as ReactElement & { ref?: React.Ref<HTMLElement> }
  ).ref;

  const wiredTrigger = cloneElement(
    trigger as ReactElement<{
      onClick?: (e: unknown) => void;
      "aria-haspopup"?: string;
      "aria-expanded"?: boolean;
      "aria-controls"?: string;
      ref?: React.Ref<HTMLElement>;
    }>,
    {
      "aria-haspopup": "dialog",
      "aria-expanded": open,
      ...(open ? { "aria-controls": panelId } : {}),
      onClick: (e: unknown) => {
        (
          trigger.props as { onClick?: (e: unknown) => void }
        ).onClick?.(e);
        setOpen((v) => !v);
      },
      ref: (node: HTMLElement | null) => {
        triggerRef.current = node;
        if (typeof existingRef === "function") existingRef(node);
        else if (existingRef && typeof existingRef === "object") {
          (existingRef as React.MutableRefObject<HTMLElement | null>).current = node;
        }
      },
    },
  );

  const onPanelKeyDown = (e: KeyboardEvent<HTMLDivElement>): void => {
    if (e.key === "Escape") {
      e.preventDefault();
      close();
    }
  };

  return (
    <span ref={wrapperRef} className={cx("relative", fullWidth ? "flex w-full" : "inline-flex")}>
      {wiredTrigger}
      {open && (
        <div
          ref={panelRef}
          id={panelId}
          role="dialog"
          aria-label={label}
          data-testid="popover"
          data-side={side}
          tabIndex={-1}
          onKeyDown={onPanelKeyDown}
          className={cx(
            "card absolute left-0 z-30 min-w-[180px] p-2 anim-pal-scale focus:outline-none",
            side === "top" ? "bottom-full mb-1" : "top-full mt-1",
            className,
          )}
        >
          {children}
        </div>
      )}
    </span>
  );
}
