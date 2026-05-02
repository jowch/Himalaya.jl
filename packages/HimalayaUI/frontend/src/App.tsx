import "./styles.css";
import { useEffect } from "react";
import { useQueryClient } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { OnboardingFlow } from "./components/OnboardingFlow";
import { useAppState } from "./state";
import { handleCurationEvent } from "./lib/sseSubscriber";

/**
 * App — root. Two concerns: the persistent workspace shell and the
 * onboarding overlay (rendered only when username is unset).
 */
export function App(): JSX.Element {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);

  useEffect(() => {
    const es = new EventSource("/api/events");
    es.addEventListener("curation", (e) => {
      handleCurationEvent((e as MessageEvent).data as string, { username, qc });
    });
    return () => es.close();
  }, [username, qc]);

  return (
    <>
      <AppShell />
      <OnboardingFlow />
    </>
  );
}
