import { useEffect } from "react";
import { Routes, Route, Navigate, useParams, useNavigate } from "react-router-dom";
import { Outlet } from "react-router-dom";
import { useAppState } from "../../state";
import { useInteraction } from "../interaction/registry";
import { core } from "../interaction/core";
import { TopNav } from "./TopNav";
import { ExperimentShell } from "./ExperimentShell";
import { ExperimentsHomePage } from "../pages/ExperimentsHomePage";
import { NewExperimentPage } from "../pages/NewExperimentPage";
import { ExperimentEntry } from "./ExperimentEntry";
import { ConfigurationPage } from "../pages/ConfigurationPage";
import { GroupingReviewPage } from "../pages/GroupingReviewPage";
import { IndexSlugRedirect } from "./IndexSlugRedirect";
import { SampleRedirect } from "./SampleRedirect";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";
import { useStateFromUrl } from "../../hooks/useStateFromUrl";
import { InteractionDock } from "../interaction/InteractionDock";
import { useKeyboardLayer } from "../interaction/useKeyboardLayer";
import { WorkspacePage } from "../../workspace/WorkspacePage";

/** Thin wrapper that reads :id from the route and passes it to GroupingReviewPage. */
function GroupingReviewRoute(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const experimentId = id ? Number(id) : 0;
  return (
    <GroupingReviewPage
      // Remount per experiment so ephemeral selection/cursor state can't carry a
      // stale multi-selection across an experiment switch.
      key={experimentId}
      experimentId={experimentId}
      onBack={() => navigate(`/experiments/${experimentId}`)}
    />
  );
}

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

/**
 * AppShell — the single unified shell wrapping every route (T3.2).
 * Renders TopNav (the unified nav bar) above the matched child.
 * Replaces the three legacy shells (deleted in T3.2).
 */
function AppShell(): JSX.Element {
  useKeyboardLayer(); // one window listener for the whole app

  // Register the shell-level global shortcuts once.  These are scope-exempt
  // (bypass the bare-key inPageScope guard) so `/`, `?`, and ⌘K work from
  // anywhere in the app, not just inside [data-interaction-scope] elements.
  useEffect(() => {
    useInteraction.getState().setShellActions([
      core("find", {
        run: () => {
          const s = useAppState.getState();
          s.openNavModal(s.activeExperimentId === undefined ? "experiment" : "sample");
        },
      }),
      core("help", { run: () => useAppState.getState().openHelpOverlay() }),
    ]);
  }, []);

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink"
    >
      <TopNav />
      <main className="flex-1 min-h-0 overflow-auto">
        <Outlet />
      </main>
      <InteractionDock />
    </div>
  );
}

export function AppRoutes(): JSX.Element {
  return (
    <Routes>
      <Route element={<AppShell />}>
        <Route path="/series/*" element={<Navigate to="/experiments" replace />} />
        <Route path="/inspect/*" element={<Navigate to="/experiments" replace />} />
        <Route path="/experiments" element={<ExperimentsHomePage />} />
        <Route path="/experiments/new" element={<NewExperimentPage />} />
        <Route path="/experiments/new/config" element={<ConfigurationPage />} />
        <Route path="/experiments/:id/grouping" element={<GroupingReviewRoute />} />
        <Route path="/experiments/:id" element={<ExperimentShell />}>
          <Route index element={<ExperimentEntry />} />
          <Route path="corpus" element={<Navigate to=".." replace />} />
          <Route path="config" element={<ConfigurationPage />} />
        </Route>
        <Route path="*" element={<PageBody />} />
      </Route>
      <Route path="/experiments/:id/workspace/:sampleId?" element={<WorkspacePage />} />
      <Route path="/sample/:sampleId/*" element={<SampleRedirect />} />
      <Route path="/" element={<Navigate to="/experiments" replace />} />
      <Route path="/index" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment/:sample" element={<IndexSlugRedirect />} />
      <Route path="/samples" element={<Navigate to="/experiments" replace />} />
      <Route path="/experiments/:eid/compare/*" element={<Navigate to="/experiments" replace />} />
      <Route path="/compare/all/*" element={<Navigate to="/experiments" replace />} />
    </Routes>
  );
}
