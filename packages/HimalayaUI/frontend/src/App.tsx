import "./styles.css";
import { useEffect, useRef } from "react";
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

  // Keep the latest username in a ref so the EventSource listener always
  // reads the current value without the EventSource being recreated on every
  // onboarding transition (undefined → set). The EventSource lifetime is
  // bound to mount/unmount only; qc is stable from QueryClientProvider.
  const usernameRef = useRef<string | undefined>(username);
  useEffect(() => {
    usernameRef.current = username;
  }, [username]);

  useEffect(() => {
    const es = new EventSource("/api/events");
    es.addEventListener("curation", (e) => {
      handleCurationEvent((e as MessageEvent).data as string, {
        username: usernameRef.current,
        qc,
      });
    });
    return () => es.close();
  }, [qc]); // qc is stable; effective deps = [] for EventSource lifetime

  return (
    <>
      <AppShell />
      <OnboardingFlow />
    </>
  );
}
