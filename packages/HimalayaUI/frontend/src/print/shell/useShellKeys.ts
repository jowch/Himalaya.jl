import { useEffect } from "react";
import { useAppState } from "../../state";
import { suppressGlobalKeys } from "../../lib/keys";

function inDialog(e: KeyboardEvent): boolean {
  const t = e.target;
  if (t instanceof Element && t.closest('[role="dialog"], [data-keys-trap]') !== null) return true;
  return document.querySelector('[role="dialog"][aria-modal="true"]') !== null;
}

export function useShellKeys(openHelp: () => void): void {
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      const mod = e.metaKey || e.ctrlKey;
      if (mod && !e.altKey && !e.shiftKey && e.key.toLowerCase() === "k") {
        if (e.defaultPrevented) return;
        const s = useAppState.getState();
        if (s.paletteOpen) {
          e.preventDefault();
          s.closePalette();
          return;
        }
        if (inDialog(e)) return;
        e.preventDefault();
        s.openPalette();
        return;
      }
      if (mod || e.altKey || suppressGlobalKeys(e)) return;
      if (e.key === "/") {
        e.preventDefault();
        useAppState.getState().openPalette();
      } else if (e.key === "?") {
        e.preventDefault();
        openHelp();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [openHelp]);
}
