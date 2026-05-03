import "./styles.css";
import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { OnboardingFlow } from "./components/OnboardingFlow";
import { handleCurationEvent } from "./lib/sseSubscriber";
import { getClientId } from "./lib/clientId";

/**
 * App — root. Two concerns: the persistent workspace shell and the
 * onboarding overlay (rendered only when username is unset).
 */
export function App(): JSX.Element {
  const qc = useQueryClient();
  const clientId = getClientId(); // stable for the tab session (sessionStorage)

  useEffect(() => {
    const es = new EventSource("/api/events");
    es.addEventListener("curation", (e) => {
      handleCurationEvent((e as MessageEvent).data as string, {
        clientId,
        qc,
      });
    });
    return () => es.close();
  }, [qc, clientId]); // both stable; effective deps = [] for EventSource lifetime

  return (
    <>
      <AppShell />
      <OnboardingFlow />
    </>
  );
}
