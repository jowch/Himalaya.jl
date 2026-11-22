import { useCallback, useState } from "react";
import { Routes, Route, Navigate, Outlet, useLocation } from "react-router-dom";
import { useAppState } from "../../state";
import { TopNav } from "./TopNav";
import { ExperimentsHomePage } from "../pages/ExperimentsHomePage";
import { NewExperimentPage } from "../pages/NewExperimentPage";
import { ExperimentEntry } from "./ExperimentEntry";
import { IndexSlugRedirect } from "./IndexSlugRedirect";
import { SampleRedirect } from "./SampleRedirect";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";
import { useStateFromUrl } from "../../hooks/useStateFromUrl";
import { useShellKeys } from "./useShellKeys";
import { ShellPalette } from "./ShellPalette";
import { HelpOverlay } from "../../workspace/ui/HelpOverlay";
import { shellHelp } from "./shellHelp";
import { WorkspacePage } from "../../workspace/WorkspacePage";

/**
 * PageBody — the `*` catch-all body. The URL→stale classifier
 * (`useStateFromUrl`) is mounted HERE, inside the catch-all element, NOT in
 * any layout shell. PageBody renders only when no other route matched, so the
 * classifier runs only on a genuinely-unmatched path — mounting it in the
 * shell would wrongly park a valid corpus route on StaleUrlPage.
 *
 * Any path that reaches `*` is, by definition, one no other route matched —
 * `useStateFromUrl` classifies it as a stale unknown_path and sets
 * `staleUrlContext` inside its effect. On the first render that effect hasn't
 * run yet (so `staleUrlContext` is still null); we render the neutral
 * ResolvingFallback for that one frame rather than navigating away.
 */
function PageBody(): JSX.Element {
  useStateFromUrl();
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  if (staleUrlContext !== null) {
    return <StaleUrlPage staleUrlContext={staleUrlContext} />;
  }
  // `resolving`, or the one pre-effect frame before useStateFromUrl sets stale.
  return <ResolvingFallback />;
}

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

export function AppRoutes(): JSX.Element {
  return (
    <Routes>
      <Route element={<AppShell />}>
        <Route path="/inspect/*" element={<Navigate to="/experiments" replace />} />
        <Route path="/experiments" element={<ExperimentsHomePage />} />
        <Route path="/experiments/new" element={<NewExperimentPage />} />
        <Route path="/experiments/new/config" element={<Navigate to="/experiments/new" replace />} />
        <Route path="/experiments/:id" element={<ExperimentEntry />} />
        {["corpus", "config", "grouping"].map((retired) => (
          <Route key={retired} path={`/experiments/:id/${retired}`} element={<ExperimentEntry />} />
        ))}
        <Route path="*" element={<PageBody />} />
      </Route>
      <Route path="/experiments/:id/workspace/:sampleId?" element={<WorkspacePage />} />
      <Route path="/sample/:sampleId/*" element={<SampleRedirect />} />
      <Route path="/" element={<Navigate to="/experiments" replace />} />
      <Route path="/index" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment/:sample" element={<IndexSlugRedirect />} />
      <Route path="/samples" element={<Navigate to="/experiments" replace />} />
    </Routes>
  );
}
