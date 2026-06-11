import { useEffect } from "react";
import { useLocation, useMatch, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { suppressGlobalKeys } from "../lib/keys";
import { useCorpusSamples } from "../queries";
import { resolveRouteSampleStatus } from "./useSyncActiveSampleFromRoute";
import { useExperimentSiblings } from "./useExperimentSiblings";
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
 *
 * F5: on the focus route the step derives its sibling list from the CORPUS
 * cache via the shared `useExperimentSiblings` derivation — the SAME list the
 * topbar stepper renders, so the two can never disagree. The
 * `samplesInExperiment` parameter is gated on `activeExperimentId`, which only
 * the NavModal picker and recoverFromStaleUrl ever write; direct-visit and
 * door-entry flows never set it, which is exactly why the shortcut used to be
 * dead on /sample/:id. Cold corpus cache or an unknown active sample degrade
 * to an empty derivation, so the shortcut no-ops gracefully there.
 */
export function useGlobalShortcuts(samplesInExperiment: Sample[] | undefined): void {
  const navigate = useNavigate();
  const { pathname } = useLocation();
  const { prev: prevSibling, next: nextSibling } = useExperimentSiblings();

  // F-STALEURL honesty gate (M1): a bogus /sample/:id never seeds the store,
  // so mid-session the previous VALID activeSampleId survives — the sibling
  // derivation above would then step relative to a sample the URL does not
  // show, out of a "Sample not found" page. Judge the route param against the
  // corpus with the SAME predicate the topbar stepper uses to hide itself
  // (resolveRouteSampleStatus); gate on "unknown" only — "pending" (cold
  // cache) already no-ops via the empty derivation.
  const sampleMatch = useMatch("/sample/:sampleId");
  const corpusQ = useCorpusSamples();
  const routeStatus = resolveRouteSampleStatus(
    sampleMatch?.params.sampleId,
    corpusQ.data,
  );

  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent): void => {
      // Typing contexts + open modals suppress. NOTE: the helper deliberately
      // does not swallow modifier chords — ⌘K below must still pass.
      if (suppressGlobalKeys(e)) return;

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
        // Focus route: the URL is the source of truth (one-way URL->store sync,
        // useSyncActiveSampleFromRoute) — a store-only write is reverted on the
        // next render. Navigate the URL through the corpus-derived sibling list
        // the topbar stepper uses (F5, shared derivation). prev/next are
        // undefined at the ends AND when the derivation is unresolved (cold
        // corpus cache, unknown STORE sample) — both no-op. The routeStatus
        // gate covers the remaining case: a bogus URL over a stale-valid store
        // sample (see the F-STALEURL note above).
        if (pathname.startsWith("/sample/")) {
          if (routeStatus === "unknown") return;
          const target = e.key === "." ? nextSibling : prevSibling;
          if (target) {
            e.preventDefault();
            navigate(`/sample/${target.id}`);
          }
          return;
        }
        // Corpus surfaces stay store-driven off the per-experiment list.
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
  }, [samplesInExperiment, navigate, pathname, prevSibling, nextSibling, routeStatus]);
}
