import { useEffect } from "react";
import { useAppState } from "../../state";
import { helpKeyFromEmptyField, inDialog, suppressGlobalKeys } from "../../lib/keys";

export function useShellKeys(openHelp: () => void): void {
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      const mod = e.metaKey || e.ctrlKey;
      if (mod && !e.altKey && !e.shiftKey && e.key.toLowerCase() === "k") {
        if (e.defaultPrevented || inDialog(e)) return;
        e.preventDefault();
        const s = useAppState.getState();
        s.openNavModal(s.activeExperimentId === undefined ? "experiment" : "sample");
        return;
      }
      if (helpKeyFromEmptyField(e)) {
        e.preventDefault();
        openHelp();
        return;
      }
      if (mod || e.altKey || suppressGlobalKeys(e)) return;
      if (e.key === "/") {
        e.preventDefault();
        const s = useAppState.getState();
        s.openNavModal(s.activeExperimentId === undefined ? "experiment" : "sample");
      } else if (e.key === "?") {
        e.preventDefault();
        openHelp();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [openHelp]);
}
