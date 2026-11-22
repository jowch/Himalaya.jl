import { useEffect } from "react";

const FOCUSABLE =
  'a[href],button:not([disabled]),input:not([disabled]),select:not([disabled]),textarea:not([disabled]),[tabindex]:not([tabindex="-1"])';

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
    if (!container.contains(document.activeElement)) focusInto();
    const onKeyDown = (e: KeyboardEvent): void => {
      if (e.key !== "Tab" || e.defaultPrevented) return;
      const current = document.activeElement;
      if (!current || !container.contains(current) || current === container) {
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