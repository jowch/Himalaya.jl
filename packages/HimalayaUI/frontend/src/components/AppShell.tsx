import { useEffect } from "react";
import { useAppState } from "../state";
import { useSamples } from "../queries";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { IndexPage } from "../pages/IndexPage";
import { ComparePage } from "../pages/ComparePage";
import { InspectPage } from "../pages/InspectPage";
import { useGlobalShortcuts } from "../hooks/useGlobalShortcuts";

/**
 * AppShell — top-level layout with the grain background, the app header,
 * and the active page body. Owns global keyboard shortcuts.
 */
export function AppShell(): JSX.Element {
  const activePage  = useAppState((s) => s.activePage);
  const theme       = useAppState((s) => s.theme);
  const experimentId = useAppState((s) => s.activeExperimentId);

  // Apply theme to <html> so our CSS can key off `html.theme-light`.
  useEffect(() => {
    const cls = theme === "light" ? "theme-light" : "";
    document.documentElement.className = cls;
    return () => { document.documentElement.className = ""; };
  }, [theme]);

  const samplesQ = useSamples(experimentId ?? 0);
  useGlobalShortcuts(experimentId === undefined ? undefined : samplesQ.data);

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full max-w-[1600px] mx-auto flex flex-col min-h-0 relative"
    >
      <AppHeader />
      {/* Page-nav rocker sits in its own row, where the per-page title
          used to live. The page title now lives in the plot card's top
          strip on the Index page. */}
      <div className="shrink-0 flex justify-center pt-1 pb-2">
        <TabRocker />
      </div>
      {activePage === "index"   && <IndexPage />}
      {activePage === "inspect" && <InspectPage />}
      {activePage === "compare" && <ComparePage />}
      <NavModal />
    </div>
  );
}
