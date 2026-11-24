import { useCallback, useState } from "react";
import { Routes, Route, Navigate, Outlet, useLocation, useParams } from "react-router";
import { TopNav } from "./TopNav";
import { ExperimentsHomePage } from "../pages/ExperimentsHomePage";
import { NewExperimentPage } from "../pages/NewExperimentPage";
import { ExperimentEntry } from "./ExperimentEntry";
import { LinkGone } from "./LinkGone";
import { useShellKeys } from "./useShellKeys";
import { ShellPalette } from "./ShellPalette";
import { HelpOverlay } from "../ui/HelpOverlay";
import { shellHelp } from "./shellHelp";
import { WorkspacePage } from "../workspace/WorkspacePage";

export const PAGE_GONE = "That page no longer exists.";

function AppShell(): JSX.Element {
  const [helpOpen, setHelpOpen] = useState(false);
  const openHelp = useCallback(() => setHelpOpen(true), []);
  useShellKeys(openHelp);
  const { pathname } = useLocation();

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink"
    >
      <TopNav />
      <main className="flex-1 min-h-0 overflow-auto">
        <Outlet />
      </main>
      <ShellPalette openHelp={openHelp} />
      {helpOpen && <HelpOverlay sections={shellHelp(pathname)} onClose={() => setHelpOpen(false)} />}
    </div>
  );
}

function WorkspaceRoute(): JSX.Element {
  const { id } = useParams();
  return <WorkspacePage key={id} />;
}

export function AppRoutes(): JSX.Element {
  return (
    <Routes>
      <Route element={<AppShell />}>
        <Route path="/experiments" element={<ExperimentsHomePage />} />
        <Route path="/experiments/new" element={<NewExperimentPage />} />
        <Route path="/experiments/:id" element={<ExperimentEntry />} />
        <Route path="*" element={<LinkGone message={PAGE_GONE} />} />
      </Route>
      <Route path="/experiments/:id/workspace/:sampleId?" element={<WorkspaceRoute />} />
      <Route path="/" element={<Navigate to="/experiments" replace />} />
    </Routes>
  );
}
