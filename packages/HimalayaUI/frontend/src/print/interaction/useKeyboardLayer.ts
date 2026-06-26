import { useEffect } from "react";
import type { Action } from "./types";
import { useInteraction } from "./registry";
import { isTyping, matchesKeys } from "./matchKey";
import { isNativeInteractiveTarget } from "../../lib/keys";

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
      // Focus in a text field → the field owns the keyboard (native undo/redo/
      // select-all, typing). No registry shortcut fires through it. (Was bare-only;
      // widened so chorded editing shortcuts like Mod+z reach native text-undo,
      // not the page's undo action.)
      if (isTyping(e.target)) return;

      const { actions, shellActions } = useInteraction.getState();
      // shellActions are scope-exempt (global shortcuts like /, ?, ⌘K).
      // page actions are scope-gated (WCAG 2.1.4: bare key fires only inside interaction scope).
      for (const [a, isShell] of [
        ...shellActions.map((a) => [a, true] as const),
        ...actions.map((a) => [a, false] as const),
      ]) {
        if (!a.keys || !matchesKeys(e, a.keys)) continue;
        // WCAG 2.1.4: bare single-key page actions fire only inside the page scope.
        // Shell actions bypass this guard (they are always global).
        if (!isShell && !inPageScope(e.target)) continue;
        // Enter on a native interactive control (button/link/input/select/textarea/
        // contenteditable) activates THAT control — don't let the page's Enter action
        // (openFocus/"Apply") hijack it. Space is already covered (it is bare → scope-gated).
        if (e.key === "Enter" && isNativeInteractiveTarget(e)) continue;
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
