import "./styles.css";
import { AppShell } from "./components/AppShell";
import { OnboardingFlow } from "./components/OnboardingFlow";

/**
 * App — root. Two concerns: the persistent workspace shell and the
 * onboarding overlay (rendered only when username is unset).
 */
export function App(): JSX.Element {
  return (
    <>
      <AppShell />
      <OnboardingFlow />
    </>
  );
}
