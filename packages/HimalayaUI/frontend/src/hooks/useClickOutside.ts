import { useEffect, useRef, type RefObject } from "react";

/**
 * While `active`, call `onOutside` when a pointerdown lands outside `ref`.
 * Listener is bound only while active and rebinds only when active/ref change
 * (the callback is held in a ref), so passing an inline `() => setOpen(false)`
 * does not re-subscribe on every render.
 */
export function useClickOutside(
  active: boolean,
  ref: RefObject<HTMLElement | null>,
  onOutside: () => void,
): void {
  const cb = useRef(onOutside);
  cb.current = onOutside;
  useEffect(() => {
    if (!active) return;
    const onPointerDown = (e: PointerEvent): void => {
      if (ref.current && !ref.current.contains(e.target as Node)) cb.current();
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, [active, ref]);
}
