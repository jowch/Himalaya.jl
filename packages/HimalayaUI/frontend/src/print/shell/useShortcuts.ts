import { useEffect, useRef } from "react";
import { matchShortcut } from "./shortcuts";
import type { ShortcutId } from "./shortcuts";
import { suppressGlobalKeys } from "../../lib/keys";

/** A binding returns nothing when it handles the key (the event is then
 *  `preventDefault`ed), or `false` to DECLINE it — leaving the event
 *  un-prevented so another listener (e.g. TracePlate's own Escape-to-disarm,
 *  the next rung of the Esc ladder) can act on it. */
export type ShortcutHandler = (e: KeyboardEvent) => void | boolean;
export type ShortcutBindings = Partial<Record<ShortcutId, ShortcutHandler>>;

/**
 * Bind window-level keyboard shortcuts from the shared registry. The single
 * implementation every surface uses instead of hand-rolling its own keydown
 * effect — so a key always means the same action and `suppressGlobalKeys`
 * (typing fields, open popovers/modals) is honored uniformly.
 *
 * `bindings` may change every render; the handler is read through a ref so the
 * listener is bound once (no churn, always the latest closures). When a binding
 * matches, the event is `preventDefault`ed; unmatched keys pass through.
 */
export function useShortcuts(bindings: ShortcutBindings, enabled: boolean = true): void {
  const ref = useRef(bindings);
  ref.current = bindings;

  useEffect(() => {
    if (!enabled) return undefined;
    function onKeyDown(e: KeyboardEvent): void {
      if (suppressGlobalKeys(e)) return;
      const id = matchShortcut(e);
      if (id === null) return;
      const handler = ref.current[id];
      if (handler === undefined) return;
      // preventDefault only if the handler actually claimed the key. A handler
      // returning `false` declined it, so the event keeps propagating (and stays
      // un-prevented) for the next listener / Esc-ladder rung.
      if (handler(e) !== false) e.preventDefault();
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [enabled]);
}
