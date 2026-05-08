import { useEffect } from "react";
import { Routes, Route, useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { useSamples } from "../queries";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { IndexPage } from "../pages/IndexPage";
import { ComparePage } from "../pages/ComparePage";
import { ComparePageEdit } from "../pages/ComparePageEdit";
import { InspectPage } from "../pages/InspectPage";
import { useGlobalShortcuts } from "../hooks/useGlobalShortcuts";

/**
 * AppShell — top-level layout with the grain background, the app header,
 * and the active page body. Owns global keyboard shortcuts.
 *
 * Routing model: the Compare page is URL-routed (so `:eid`/`:id` survive
 * reloads — see Plan §Phase 4). Index and Inspect remain Zustand-driven
 * for now (existing behaviour). The two systems coexist via:
 *   - URL `/experiments/:eid/compare*` or `/compare/all` → render Compare
 *   - Anything else → render IndexPage / InspectPage based on `activePage`
 *
 * `TabRocker` syncs the two: clicking "Compare" navigates to the URL,
 * clicking "Index"/"Inspect" navigates back to "/" and updates `activePage`.
 */
function ZustandShellPage(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  return (
    <>
      {activePage === "index"   && <IndexPage />}
      {activePage === "inspect" && <InspectPage />}
    </>
  );
}

export function AppShell(): JSX.Element {
  const theme        = useAppState((s) => s.theme);
  const experimentId = useAppState((s) => s.activeExperimentId);
  const setActivePage = useAppState((s) => s.setActivePage);
  const location = useLocation();
  const navigate = useNavigate();

  // Apply theme to <html> so our CSS can key off `html.theme-light`.
  useEffect(() => {
    const cls = theme === "light" ? "theme-light" : "";
    document.documentElement.className = cls;
    return () => { document.documentElement.className = ""; };
  }, [theme]);

  // Sync URL → Zustand activePage. When the URL is /compare* or
  // /experiments/:eid/compare*, mark the page tab as "compare". When the
  // URL is "/", we leave activePage alone (it's already index/inspect).
  const onComparePath = location.pathname.startsWith("/compare")
    || /^\/experiments\/\d+\/compare(\/|$)/.test(location.pathname);
  useEffect(() => {
    if (onComparePath) setActivePage("compare");
  }, [onComparePath, setActivePage]);

  const activePage = useAppState((s) => s.activePage);

  // Symmetric: when activePage is "compare" but URL isn't on a compare path,
  // navigate so the URL-routed Compare pages mount. Without this, a reload
  // at "/" with persisted activePage='compare' renders the rocker but no
  // page body (issue #77).
  useEffect(() => {
    if (activePage !== "compare") return;
    if (onComparePath) return;
    const url = experimentId !== undefined
      ? `/experiments/${experimentId}/compare`
      : "/compare/all";
    navigate(url, { replace: true });
  }, [activePage, onComparePath, experimentId, navigate]);

  // When the user's activePage flips from compare back to index/inspect via
  // TabRocker, we navigate back to "/" so the Zustand shell takes over.
  // This is handled in TabRocker itself — see TabRocker.tsx.

  const samplesQ = useSamples(experimentId ?? 0);
  useGlobalShortcuts(experimentId === undefined ? undefined : samplesQ.data);

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full max-w-[1600px] mx-auto flex flex-col min-h-0 relative"
      // --chrome-h: AppHeader (h-11 = 44px) + TabRocker row (pt-1 + ~28px + pb-2 ≈ 40px).
      // Pages use this to cap grid height without repeating the sum.
      style={{ "--chrome-h": "84px" } as React.CSSProperties}
    >
      <AppHeader />
      {/* Page-nav rocker sits in its own row, where the per-page title
          used to live. The page title now lives in the plot card's top
          strip on the Index page. */}
      <div className="shrink-0 flex justify-center pt-1 pb-2">
        <TabRocker
          experimentId={experimentId}
          onNavigateAway={(target) => {
            // Leaving Compare → return to "/" so the Zustand shell renders
            // the chosen page.
            if (onComparePath && target !== "compare") navigate("/");
          }}
        />
      </div>
      <Routes>
        <Route path="/experiments/:eid/compare" element={<ComparePage />} />
        <Route path="/experiments/:eid/compare/new" element={<ComparePageEdit />} />
        <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
        <Route path="/experiments/:eid/compare/:id/edit" element={<ComparePageEdit />} />
        <Route path="/compare/all" element={<ComparePage />} />
        {/*
          Global (experiment-less) deep-link routes — mirror the experiment-
          scoped review/edit/new triple so picking a comparison from
          /compare/all lands on its review page (not the empty list). New
          drafts created from /compare/all/new have no `experiment_id`
          association on the backend (comparisons aren't FK'd to experiments
          per spec); the global picker context applies.
        */}
        <Route path="/compare/all/new" element={<ComparePageEdit />} />
        <Route path="/compare/all/:id" element={<ComparePage />} />
        <Route path="/compare/all/:id/edit" element={<ComparePageEdit />} />
        <Route path="*" element={<ZustandShellPage />} />
      </Routes>
      <NavModal />
    </div>
  );
}
