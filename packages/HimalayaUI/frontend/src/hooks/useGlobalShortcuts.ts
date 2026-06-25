import { useEffect } from "react";
import { useAppState } from "../state";
import { suppressGlobalKeys } from "../lib/keys";

/**
 * useGlobalShortcuts — the app-wide keyboard binding for the ONE gesture that
 * has no home surface: find / jump.
 *
 *   `/`, `⌘K`  — open the nav modal (experiment step if none selected, else sample)
 *
 * Everything else is surface-owned through the shortcut library
 * (`print/shell/shortcuts.ts` + `useShortcuts`), so a key means the same action
 * everywhere it appears:
 *   - the sample step `[`/`]` lives on the surfaces that have it (Loupe, Focus),
 *     NOT here. KEYS-LIB step 2 retired the old global `,`/`.` stepper — `[`/`]`
 *     is the single sample-step gesture now, and the contact sheet navigation is
 *     owned at the page level.
 *   - the exposure/candidate arrows, drop/keep/representative, undo/redo and
 *     reorder gestures are all bound by their owning page via `useShortcuts`.
 *
 * R0a (#221): the `T` theme-toggle shortcut is retired with the dark theme —
 * "The Print" is the single identity, so there is no theme to toggle.
 *
 * The find chord is suppressed when typing in an input/textarea (via
 * `suppressGlobalKeys`); the modal owns its own Esc/Enter/Backspace. NOTE: the
 * suppression helper deliberately does NOT swallow modifier chords, so ⌘K still
 * passes through from a non-editing target.
 */
export function useGlobalShortcuts(): void {
  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent): void => {
      if (suppressGlobalKeys(e)) return;

      // `/` or `⌘K` → nav modal (the find/jump gesture).
      const isModK = (e.key === "k" || e.key === "K") && (e.metaKey || e.ctrlKey);
      const isSlash = e.key === "/" && !e.metaKey && !e.ctrlKey && !e.altKey;
      if (isModK || isSlash) {
        e.preventDefault();
        const s = useAppState.getState();
        const step = s.activeExperimentId === undefined ? "experiment" : "sample";
        s.openNavModal(step);
        return;
      }

      // `?` → keyboard shortcut overlay (helpOverlay). Layout-robust: `?` is
      // Shift+/ on US keyboards but `eventCombo` emits a stable `"?"` token for
      // e.key === "?" regardless of layout (mirrors the shortcuts.ts special case).
      if (e.key === "?") {
        e.preventDefault();
        useAppState.getState().openHelpOverlay();
      }
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, []);
}
