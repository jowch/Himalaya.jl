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

/** An open dialog/popover (`[role="dialog"]` — ModalShell, Popover) or inline-edit
 *  trap (`[data-keys-trap]`) OWNS the keyboard while focus is inside it. Because
 *  those overlays focus-trap, the event target sits within them — so this is
 *  target-scoped and RACE-FREE (no global `[aria-modal]` querySelector, which
 *  loses to the dialog-close microtask between the dialog's and the window's
 *  listeners). Returns true → the keyboard layer must claim nothing. */
function inOverlay(target: EventTarget | null): boolean {
  return target instanceof HTMLElement && target.closest('[role="dialog"],[data-keys-trap]') !== null;
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
      // A focus-trapped dialog/popover/inline-edit overlay owns the keyboard —
      // claim nothing through it. Critically this stops scope-exempt arrows from
      // navigating the surface BEHIND an open modal (ModalShell preventDefaults
      // only Escape, never arrows).
      if (inOverlay(e.target)) return;

      const { actions, shellActions, arrowHandler } = useInteraction.getState();

      // Arrow navigation drives the active surface GLOBALLY — scope-exempt,
      // because arrows are not WCAG-2.1.4 character shortcuts, so they need no
      // focus container. This is what stops the page from scrolling when focus
      // sits outside the surface (on body, a dock button, blank space). Widgets
      // that legitimately consume arrows (SegmentedControl, listboxes) call
      // preventDefault and were already returned at the top; text fields are
      // caught by isTyping. The page claims the arrows it handles via
      // preventDefault; an unclaimed arrow (e.g. Alt+Arrow reorder) falls through
      // to the scope-gated action loop below.
      if (arrowHandler && e.key.startsWith("Arrow")) {
        arrowHandler(e);
        if (e.defaultPrevented) return;
      }
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
        if (e.key === "Enter" && !e.metaKey && !e.ctrlKey && !e.altKey && isNativeInteractiveTarget(e)) continue;
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
