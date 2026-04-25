import { useEffect } from "react";
import { useAppState } from "../state";
import type { Sample } from "../api";

/**
 * useGlobalShortcuts — wires the keyboard shortcuts described in the plan.
 *
 *   `/`, `⌘K`  — open the nav modal (experiment step if none selected, else sample)
 *   `,`, `.`   — previous / next sample within the current experiment
 *   `T`        — toggle theme
 *
 * Shortcuts are suppressed when typing in an input/textarea. The modal itself
 * owns its own Esc/Enter/Backspace behavior; we don't touch them here.
 */
export function useGlobalShortcuts(samplesInExperiment: Sample[] | undefined): void {
  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent): void => {
      const t = e.target as HTMLElement | null;
      const editing = t && (
        t.tagName === "INPUT" || t.tagName === "TEXTAREA" ||
        (t as HTMLElement).isContentEditable
      );
      if (editing) return;

      // `/` or `⌘K` → nav modal
      if ((e.key === "k" || e.key === "K") && (e.metaKey || e.ctrlKey)) {
        e.preventDefault();
        const s = useAppState.getState();
        const step = s.activeExperimentId === undefined ? "experiment" : "sample";
        s.openNavModal(step);
        return;
      }
      if (e.key === "/" && !e.metaKey && !e.ctrlKey && !e.altKey) {
        e.preventDefault();
        const s = useAppState.getState();
        const step = s.activeExperimentId === undefined ? "experiment" : "sample";
        s.openNavModal(step);
        return;
      }

      // T → theme toggle
      if ((e.key === "t" || e.key === "T") && !e.metaKey && !e.ctrlKey && !e.altKey) {
        const s = useAppState.getState();
        s.setTheme(s.theme === "dark" ? "light" : "dark");
        return;
      }

      // `,` / `.` → prev / next sample within experiment (no wrap)
      if ((e.key === "," || e.key === ".") && !e.metaKey && !e.ctrlKey && !e.altKey) {
        const samples = samplesInExperiment ?? [];
        if (samples.length === 0) return;
        const cur = useAppState.getState().activeSampleId;
        const idx = cur === undefined ? -1 : samples.findIndex((s) => s.id === cur);
        const step = e.key === "." ? +1 : -1;
        const nextIdx = Math.max(0, Math.min(samples.length - 1, (idx === -1 ? 0 : idx + step)));
        const next = samples[nextIdx];
        if (next && next.id !== cur) {
          e.preventDefault();
          useAppState.getState().setActiveSample(next.id);
        }
      }
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [samplesInExperiment]);
}
