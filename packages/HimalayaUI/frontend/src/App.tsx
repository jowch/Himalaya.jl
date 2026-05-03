import "./styles.css";
import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { OnboardingFlow } from "./components/OnboardingFlow";
import { ToastContainer } from "./components/ui/Toast";
import { InfrastructureBanner } from "./components/InfrastructureBanner";
import { handleRemoteEvent } from "./lib/queue/replayCoordinator";
import { attachPersistence } from "./lib/queue/persistence";
import type { SseEvent } from "./lib/queue/types";

/**
 * App — root. Two concerns: the persistent workspace shell and the
 * onboarding overlay (rendered only when username is unset).
 */
export function App(): JSX.Element {
  const qc = useQueryClient();
  const mc = qc.getMutationCache();

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

  // Mirror pending mutation queue to sessionStorage so a tab reload could
  // rehydrate it (rehydrate wiring deferred — see M3 follow-up). attach is
  // safe to mount unconditionally; subscription is cheap.
  useEffect(() => {
    return attachPersistence(mc);
  }, [mc]);

  return (
    <>
      <AppShell />
      <OnboardingFlow />
      <ToastContainer />
      <InfrastructureBanner />
    </>
  );
}
