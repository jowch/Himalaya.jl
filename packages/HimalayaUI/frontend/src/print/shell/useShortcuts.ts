import { useEffect, useRef } from "react";
import { matchShortcut } from "./shortcuts";
import type { ShortcutId } from "./shortcuts";
import { suppressGlobalKeys } from "../../lib/keys";

export type ShortcutBindings = Partial<Record<ShortcutId, (e: KeyboardEvent) => void>>;

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
      e.preventDefault();
      handler(e);
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [enabled]);
}
