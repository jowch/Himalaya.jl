import { useEffect, useRef, type MouseEvent } from "react";

export function releaseMouseFocus(e: MouseEvent<HTMLElement>): void {
  if (e.detail > 0) e.currentTarget.blur();
}

export function useRovingFocus<T extends HTMLElement>(currentKey: unknown) {
  const ref = useRef<T | null>(null);
  useEffect(() => {
    const root = ref.current;
    const active = document.activeElement;
    if (!root || !(active instanceof HTMLElement) || !root.contains(active)) return;
    const slot = active.dataset.roving;
    if (slot === undefined || active.tabIndex === 0) return;
    root.querySelector<HTMLElement>(`[data-roving="${slot}"][tabindex="0"]`)?.focus({ preventScroll: true });
  }, [currentKey]);
  return ref;
}
