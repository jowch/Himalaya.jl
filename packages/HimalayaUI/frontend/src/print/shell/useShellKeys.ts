import { useEffect } from "react";
import { useAppState } from "../../state";
import { suppressGlobalKeys } from "../../lib/keys";

export function useShellKeys(openHelp: () => void): void {
  useEffect(() => {
    const onKey = (e: KeyboardEvent): void => {
      if (suppressGlobalKeys(e)) return;
      const mod = e.metaKey || e.ctrlKey;
      const find = (mod && !e.altKey && !e.shiftKey && e.key.toLowerCase() === "k") || (!mod && !e.altKey && e.key === "/");
      if (find) {
        e.preventDefault();
        const s = useAppState.getState();
        s.openNavModal(s.activeExperimentId === undefined ? "experiment" : "sample");
      } else if (!mod && !e.altKey && e.key === "?") {
        e.preventDefault();
        openHelp();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [openHelp]);
}
