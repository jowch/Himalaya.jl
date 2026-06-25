import { useEffect } from "react";
import type { Action } from "./types";
import { useInteraction } from "./registry";
import { isBareKey, isTyping, matchesKeys } from "./matchKey";

function inPageScope(target: EventTarget | null): boolean {
  if (!(target instanceof HTMLElement)) return true; // window/body-level keydown
  if (target === document.body) return true;
  return target.closest("[data-interaction-scope],[data-cursored]") !== null;
}

function isEnabled(a: Action): boolean {
  return a.enabled ? a.enabled() : true;
}

/** The single rung-3 keyboard entry point. Mounted once in the app shell. */
export function useKeyboardLayer(): void {
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      // Rung 1/2 already handled it (the microtask-race fix — guard the signal,
      // never querySelector('[aria-modal]')).
      if (e.defaultPrevented) return;
      if (isTyping(e.target) && isBareKey(e)) return;

      const { actions } = useInteraction.getState();
      for (const a of actions) {
        if (!a.keys || !matchesKeys(e, a.keys)) continue;
        // WCAG 2.1.4: bare single-key actions fire only inside the page scope.
        if (isBareKey(e) && !inPageScope(e.target)) continue;
        if (!isEnabled(a)) return; // matched but inert — claim nothing, swallow nothing
        e.preventDefault();
        a.run(e);
        return;
      }
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, []);
}
