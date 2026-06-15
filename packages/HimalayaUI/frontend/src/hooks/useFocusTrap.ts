import { useEffect } from "react";

const FOCUSABLE =
  'a[href],button:not([disabled]),input:not([disabled]),select:not([disabled]),textarea:not([disabled]),[tabindex]:not([tabindex="-1"])';

/**
 * Traps keyboard focus within `containerRef` while `active` is true
 * (WAI-ARIA APG modal-dialog pattern).
 *
 * - On activation, moves focus INTO the container — to its first focusable
 *   descendant, or to the container itself (given tabIndex=-1) if none —
 *   unless the consumer already placed focus inside. React `autoFocus` fires
 *   during commit (before this effect) and is respected; a PARENT's own focus
 *   effect (NavModal's input focus) instead runs AFTER this child-mounted
 *   trap and simply overwrites its initial placement — both end up inside.
 * - One trap at a time: two simultaneously-active traps would ping-pong via
 *   their focusin guards. Unreachable today (suppressGlobalKeys blocks modal
 *   openers while a dialog is up); add stacking support before nesting modals.
 * - Intercepts Tab/Shift+Tab at the DOCUMENT level, so the trap holds even
 *   when focus has escaped the container: boundary presses wrap, and a press
 *   with focus outside is driven back to the first (last for Shift) focusable.
 *   Focusables are recomputed per keydown — modal contents change. A Tab the
 *   consumer already `preventDefault`ed is left alone (NavModal Tab-commit).
 * - A document `focusin` guard pulls focus back in if it lands outside.
 * - Restores focus to the previously-focused element on deactivation.
 *   Effect setup/cleanup are symmetric, so StrictMode double-activation is
 *   safe (prevFocus is re-captured per activation).
 */
export function useFocusTrap(
  containerRef: React.RefObject<HTMLElement | null>,
  active: boolean,
): void {
  useEffect(() => {
    if (!active) return;
    const container = containerRef.current;
    if (!container) return;

    const prevFocus = document.activeElement as HTMLElement | null;

    const focusables = (): HTMLElement[] =>
      Array.from(container.querySelectorAll<HTMLElement>(FOCUSABLE));

    const focusInto = (fromEnd = false): void => {
      const list = focusables();
      const target = fromEnd ? list[list.length - 1] : list[0];
      if (target) { target.focus(); return; }
      if (!container.hasAttribute("tabindex")) container.tabIndex = -1;
      container.focus();
    };

    // Engage: move focus into the dialog unless a consumer already did.
    if (!container.contains(document.activeElement)) focusInto();

    const onKeyDown = (e: KeyboardEvent): void => {
      if (e.key !== "Tab" || e.defaultPrevented) return;
      const current = document.activeElement;
      if (!current || !container.contains(current) || current === container) {
        // Focus escaped (or sits on the container shell) — drive it back in.
        e.preventDefault();
        focusInto(e.shiftKey);
        return;
      }
      const list = focusables();
      if (list.length === 0) {
        e.preventDefault();
        focusInto();
        return;
      }
      const first = list[0]!;
      const last  = list[list.length - 1]!;
      if (e.shiftKey) {
        if (current === first) { e.preventDefault(); last.focus(); }
      } else {
        if (current === last)  { e.preventDefault(); first.focus(); }
      }
    };

    const onFocusIn = (e: FocusEvent): void => {
      if (container.contains(e.target as Node)) return;
      focusInto();
    };

    document.addEventListener("keydown", onKeyDown);
    document.addEventListener("focusin", onFocusIn);
    return () => {
      document.removeEventListener("keydown", onKeyDown);
      document.removeEventListener("focusin", onFocusIn);
      prevFocus?.focus();
    };
  }, [containerRef, active]);
}
