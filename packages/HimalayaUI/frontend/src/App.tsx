import "./styles.css";
import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { OnboardingFlow } from "./components/OnboardingFlow";
import { ToastContainer } from "./components/ui/Toast";
import { InfrastructureBanner } from "./components/InfrastructureBanner";
import { handleRemoteEvent } from "./lib/queue/replayCoordinator";
import { attachPersistence, rehydrate } from "./lib/queue/persistence";
import { resolveMutator } from "./lib/queue/mutatorRegistry";
import { exposeTestHelpers } from "./lib/queue/testHelpers";
import { showToast } from "./lib/toast";
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
      <AppShell />
      <OnboardingFlow />
      <ToastContainer />
      <InfrastructureBanner />
    </>
  );
}
