import "./styles.css";
import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppRoutes } from "./components/AppRoutes";
import { OnboardingFlow } from "./components/OnboardingFlow";
import { NavModal } from "./components/NavModal";
import { ToastContainer } from "./components/ui/Toast";
import { InfrastructureBanner } from "./components/InfrastructureBanner";
import { SeriesCommitConflictModal } from "./components/SeriesCommitConflictModal";
import { handleRemoteEvent } from "./lib/queue/replayCoordinator";
import { attachPersistence, rehydrate } from "./lib/queue/persistence";
import { attachConflictBridge } from "./lib/queue/conflictBridge";
import { resolveMutator } from "./lib/queue/mutatorRegistry";
import { exposeTestHelpers } from "./lib/queue/testHelpers";
import { showToast } from "./lib/toast";
import { useAppState } from "./state";
import type { SseEvent } from "./lib/queue/types";

/**
 * App — root. Two concerns: the persistent workspace shell and the
 * onboarding overlay (rendered only when username is unset).
 */
export function App(): JSX.Element {
  const qc = useQueryClient();
  const mc = qc.getMutationCache();

  // Expose minimal test helpers on `window.__himalayaTest` in DEV only.
  // Production bundles tree-shake this out (Vite + DEV gate).
  useEffect(() => {
    exposeTestHelpers(qc, mc);
  }, [qc, mc]);

  useEffect(() => {
    const es = new EventSource("/api/events");
    es.addEventListener("curation", (e) => {
      try {
        const parsed = JSON.parse((e as MessageEvent).data as string) as SseEvent;
        handleRemoteEvent(parsed, qc, mc);
      } catch {
        // malformed frame, ignore
      }
    });
    return () => es.close();
  }, [qc, mc]); // both stable; effective deps = [] for EventSource lifetime

  // Mirror pending mutation queue to sessionStorage so a tab reload can
  // rehydrate it. attach is safe to mount unconditionally; subscription is
  // cheap.
  useEffect(() => {
    return attachPersistence(mc);
  }, [mc]);

  // Single-source-of-truth bridge: a ConflictError on a `series_commit`
  // mutation → Zustand `pendingConflict` (the slot the
  // `SeriesCommitConflictModal` reads). Mounted once at App startup so the
  // mount sites can't race on the slot. Module-scoped last-seen tracking keeps
  // remount/HMR from re-popping the modal on a stale terminal-error mutation
  // still in the cache. See `lib/queue/conflictBridge.ts`.
  //
  // I3.6 (#177): Compare is retired, so the `comparison_save` arm of the
  // bridge is gone; only `series_commit` remains. The bridge + slot are KEPT
  // (series uses them).
  useEffect(() => {
    const setPendingConflict = useAppState.getState().setPendingConflict;
    return attachConflictBridge(mc, setPendingConflict);
  }, [mc]);

  // On mount, replay any persisted ops left over from a previous tab
  // session through their matching mutators. Server-side request-level
  // idempotency (X-Client-Op-Id) makes this safe even if the original op
  // already landed. Surfacing dropped/failed counts as toasts so the user
  // knows when edits couldn't be restored.
  useEffect(() => {
    void rehydrate(qc, resolveMutator).then(({ dropped, failed }) => {
      if (dropped > 0) {
        showToast(
          `${dropped} edits from a previous session couldn't be restored`,
          "warning",
        );
      }
      if (failed > 0) {
        showToast(
          `${failed} edits failed to replay; please retry`,
          "error",
        );
      }
    });
    // Intentionally empty deps — rehydrate runs once at mount; qc is stable.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  return (
    <>
      <AppRoutes />
      <OnboardingFlow />
      {/* I5.1 (#182): NavModal was mounted in the retired AppShell. It must
          stay mounted app-wide so the `/` + ⌘K shortcuts and StaleUrlPage
          recovery can open it from any surface. */}
      <NavModal />
      <SeriesCommitConflictModal />
      <ToastContainer />
      <InfrastructureBanner />
    </>
  );
}
