import { useEffect } from "react";

const FOCUSABLE =
  'button:not([disabled]),input:not([disabled]),select:not([disabled]),[tabindex]:not([tabindex="-1"])';

/**
 * Traps keyboard focus within `containerRef` while `active` is true.
 * Restores focus to the previously focused element on deactivation.
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

    const handler = (e: KeyboardEvent): void => {
      if (e.key !== "Tab") return;
      const focusable = Array.from(container.querySelectorAll<HTMLElement>(FOCUSABLE));
      if (focusable.length === 0) return;
      const first = focusable[0]!;
      const last  = focusable[focusable.length - 1]!;
      if (e.shiftKey) {
        if (document.activeElement === first) { e.preventDefault(); last.focus(); }
      } else {
        if (document.activeElement === last)  { e.preventDefault(); first.focus(); }
      }
    };

    container.addEventListener("keydown", handler);
    return () => {
      container.removeEventListener("keydown", handler);
      prevFocus?.focus();
    };
  }, [containerRef, active]);
}
