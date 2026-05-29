import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import type { Sample } from "../api";

/**
 * useGlobalShortcuts — wires the keyboard shortcuts described in the plan.
 *
 *   `/`, `⌘K`  — open the nav modal (experiment step if none selected, else sample)
 *   `,`, `.`   — previous / next sample within the current experiment
 *
 * R0a (#221): the `T` theme-toggle shortcut is retired with the dark theme —
 * "The Print" is the single identity, so there is no theme to toggle.
 *
 * Shortcuts are suppressed when typing in an input/textarea. The modal itself
 * owns its own Esc/Enter/Backspace behavior; we don't touch them here.
 *
 * I5.1 (#182): the ArrowLeft/Right "page-tab step" is gone with the dual-nav
 * `activePage` model — there are no legacy page tabs left to step between, so
 * the arrow keys are no longer bound here (corpus surfaces like the loupe own
 * them outright).
 *
 * The `,`/`.` step is route-aware: on the focus route (`/sample/:id`) it
 * navigates the URL (the source of truth there), elsewhere it sets the store.
 * Hence `useLocation`/`useNavigate` are read here again.
 */
export function useGlobalShortcuts(samplesInExperiment: Sample[] | undefined): void {
  const navigate = useNavigate();
  const { pathname } = useLocation();
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
          // On the focus route the URL is the source of truth (one-way
          // URL->store sync, useSyncActiveSampleFromRoute) — a store-only write
          // is reverted on the next render, which is why this shortcut read as
          // dead there. Navigate the URL instead, matching the topbar stepper.
          // Corpus surfaces stay store-driven.
          if (pathname.startsWith("/sample/")) {
            navigate(`/sample/${next.id}`);
          } else {
            useAppState.getState().setActiveSample(next.id);
          }
        }
      }
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [samplesInExperiment, navigate, pathname]);
}
