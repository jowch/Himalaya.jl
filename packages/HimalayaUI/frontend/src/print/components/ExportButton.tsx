import { useEffect, useId, useLayoutEffect, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Button, IconButton, Menu } from "../ui";
import { cx } from "../../lib/cx";


export interface ExportButtonProps {
  onCopy: () => void;
  /** Download handlers. Wiring contract: each MUST flip `pending` true
   *  synchronously (useFigureExport does) — the parked-focus repair consumes
   *  its own-activation flag at the next pending true→false edge, so a
   *  handler that never starts a pending cycle would leave the flag armed
   *  for a later (possibly foreign) cycle to consume. */
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  /** Copy unavailable (no clipboard / no PNG renderer / a render in flight). */
  copyDisabled?: boolean;
  /** PNG download unavailable (no PNG renderer / a render in flight). */
  pngDisabled?: boolean;
  /** A render is in flight (blocks SVG too). */
  pending?: boolean;
  /** Page-level gate (data not ready). ORs into every action. */
  disabled?: boolean;
  /** Why the button is disabled — when present (and `disabled`), a quiet caption
   *  renders beside the group and is wired to both actions via aria-describedby,
   *  so a disabled Export states its reason (no-data vs loading vs error) rather
   *  than going silently dead. Ignored when not disabled. */
  disabledReason?: string;
  /** Fills aria-labels: "Copy {ariaContext} to clipboard". */
  ariaContext: string;
  /** PLACEMENT ONLY. */
  className?: string;
}

/**
 * The figure-export split button (mockup: series-builder rail-foot "Copy as
 * PNG"): a bordered group with a primary **Copy** action and a `▾` chevron that
 * opens a two-item download menu (PNG / SVG). Presentational — every side
 * effect is a prop (wire it with `useFigureExport`). Owns its menu-open state,
 * outside-pointerdown dismissal, and the APG menu-button trigger half (the
 * `Field` precedent; `Menu` owns in-menu Escape + arrow-nav and moves focus
 * into the menu on open):
 * - ArrowDown on the closed trigger opens + focuses the first enabled item;
 *   ArrowUp opens + focuses the LAST (both preventDefault — no page scroll).
 * - Escape on the trigger while open closes (focus already on the trigger).
 * - Keyboard-reachable closes (Escape in the menu, item select) RETURN focus
 *   to the trigger; outside-pointerdown close does NOT — the pointer user
 *   clicked elsewhere on purpose, so we never yank focus back.
 */
export function ExportButton({
  onCopy,
  onDownloadPng,
  onDownloadSvg,
  copyDisabled = false,
  pngDisabled = false,
  pending = false,
  disabled = false,
  disabledReason,
  ariaContext,
  className,
}: ExportButtonProps): JSX.Element {
  const reasonId = useId();
  const showReason = disabled && !!disabledReason;
  const [open, setOpen] = useState(false);
  // Where Menu's mount effect puts focus on open: "first" (click/ArrowDown)
  // or "last" (ArrowUp). APG menu-button.
  const [menuFocus, setMenuFocus] = useState<"first" | "last">("first");
  const wrapRef = useRef<HTMLSpanElement | null>(null);
  const triggerRef = useRef<HTMLButtonElement | null>(null);

  // Outside-pointerdown closes the menu (mirrors Field/Popover). Bound only
  // while open.
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

  const copyOff = disabled || copyDisabled;
  const pngOff = disabled || pngDisabled;
  const svgOff = disabled || pending;

  // Close + return focus to the trigger (keyboard close paths only —
  // outside-pointerdown deliberately bypasses this; see the doc comment).
  const closeAndRefocus = (): void => {
    setOpen(false);
    triggerRef.current?.focus();
  };

  // Parked focus across the pending window (RepresentativeBox precedent):
  // selecting a download flips `pending` synchronously, which disables the
  // trigger (`disabled || pending`) — and real browsers BLUR a focused element
  // the moment it disables, so closeAndRefocus()'s focus is transient and the
  // keyboard user is dumped on <body> for the render window. Park the intent
  // in a ref set ONLY by our own menu select; a layout effect restores focus
  // to the trigger when pending transitions back to false (before paint).
  // Foreign pending flips (an export initiated elsewhere) never set the flag,
  // so they never yank focus. The flag is consumed at every pending end.
  const refocusAfterPending = useRef(false);
  const prevPending = useRef(pending);
  useLayoutEffect(() => {
    const was = prevPending.current;
    prevPending.current = pending;
    if (was && !pending) {
      if (refocusAfterPending.current) triggerRef.current?.focus();
      refocusAfterPending.current = false; // consume — lives one cycle
    }
  }, [pending]);

  // APG menu-button trigger keys: arrows open with a focus target; Escape
  // closes an open menu whose focus is still on the trigger (Menu only hears
  // Escape once focus is inside it). Enter/Space open via the native click.
  const onTriggerKeyDown = (e: KeyboardEvent<HTMLButtonElement>): void => {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault(); // arrows must not scroll the page
      setMenuFocus(e.key === "ArrowDown" ? "first" : "last");
      setOpen(true);
    } else if (e.key === "Escape" && open) {
      e.preventDefault();
      e.stopPropagation(); // innermost popup only — don't double-close a host modal
      setOpen(false); // focus is already on the trigger
    }
  };

  return (
    <span
      ref={wrapRef}
      data-testid="export-button"
      // UI-MENUTAB — APG menu-button: Tab (Shift+Tab) closes the menu and lets
      // focus proceed. Mirror the outside-pointerdown close on the FOCUS axis:
      // when focus leaves the whole widget (relatedTarget outside wrapRef, or
      // null), close. Focus moving to the trigger or between items stays inside
      // wrapRef, so a toggle-click never closes-then-reopens. Plain setOpen
      // (no refocus) — the user is already leaving.
      onBlur={(e) => {
        if (wrapRef.current && !wrapRef.current.contains(e.relatedTarget as Node | null)) {
          setOpen(false);
        }
      }}
      className={cx("relative inline-flex items-center gap-2", className)}
    >
      {showReason && (
        <span
          id={reasonId}
          role="note"
          data-testid="export-disabled-reason"
          className="text-caption text-ink-soft"
        >
          {disabledReason}
        </span>
      )}
      <span className="inline-flex items-stretch border border-hair-strong rounded overflow-hidden">
        <Button
          variant="ghost"
          data-testid="export-copy"
          aria-label={`Copy ${ariaContext} to clipboard`}
          disabled={copyOff}
          {...(showReason ? { "aria-describedby": reasonId } : {})}
          onClick={onCopy}
        >
          Copy
        </Button>
        <span className="w-px bg-hair-strong" aria-hidden="true" />
        <IconButton
          ref={triggerRef}
          label="Download formats"
          data-testid="export-menu-trigger"
          aria-haspopup="menu"
          aria-expanded={open}
          disabled={disabled || pending}
          {...(showReason ? { "aria-describedby": reasonId } : {})}
          onClick={() => {
            setMenuFocus("first"); // click-open focuses the first enabled item
            setOpen((o) => !o);
          }}
          onKeyDown={onTriggerKeyDown}
        >
          ▾
        </IconButton>
      </span>
      <Menu<"png" | "svg">
        open={open}
        options={[
          { value: "png", label: "Download as PNG", disabled: pngOff },
          { value: "svg", label: "Download as SVG", disabled: svgOff },
        ]}
        onSelect={(v) => {
          refocusAfterPending.current = true; // our activation — park focus intent
          if (v === "png") onDownloadPng();
          else onDownloadSvg();
        }}
        onClose={closeAndRefocus}
        aria-label="Download formats"
        initialFocus={menuFocus}
        className="right-0 top-full"
      />
    </span>
  );
}
