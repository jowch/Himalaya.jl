import { useEffect } from "react";
import { Routes, Route, Navigate, useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { useSamples } from "../queries";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { IndexPage } from "../pages/IndexPage";
import { ComparePage } from "../pages/ComparePage";
import { ComparePageEdit } from "../pages/ComparePageEdit";
import { InspectPage } from "../pages/InspectPage";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";
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
function PageBody(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const resolving = useAppState((s) => s.resolving);
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  if (resolving) return <ResolvingFallback />;
  if (staleUrlContext !== null) {
    return <StaleUrlPage staleUrlContext={staleUrlContext} />;
  }
  if (activePage === "index")   return <IndexPage />;
  if (activePage === "inspect") return <InspectPage />;
  // activePage === "compare" never reaches here because compare URLs are
  // matched by their explicit <Route> entries above.
  return <></>;
}

/**
 * EditToBareRedirect — Compare UX Phase B: the `/edit` URL segment is gone.
 * `/compare/:id` is the only route shape; edit is now an inline gesture, not
 * a separate URL. Old `/edit` deep-links still resolve by redirecting to the
 * bare path.
 *
 * Reads `useLocation()` (the router's internal store), NOT
 * `window.location` — under `MemoryRouter` (used in tests) JSDOM's
 * `window.location.pathname` stays at "/" and the redirect would land in the
 * wrong place. `useLocation` is correct in both browser and test envs.
 */
function EditToBareRedirect(): JSX.Element {
  const loc = useLocation();
  const here = loc.pathname.replace(/\/edit\/?$/, "");
  return <Navigate to={here} replace />;
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
  //
  // Trade-off: browser-Back from /compare/<id> to "/" bounces immediately
  // back to a compare URL (Back doesn't clear `activePage`). This is the
  // intentional fallout of the existing localStorage-persisted activePage
  // model — once the user is on Compare, "/" with no URL trail is treated
  // as a redirect target, not a meaningful destination. Users who want
  // Index/Inspect should click the rocker (which clears activePage). The
  // alternative (no redirect) was the empty-body bug we're fixing.
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
        {/* Phase B: `/edit` is gone — redirect old deep-links to the bare URL. */}
        <Route path="/experiments/:eid/compare/:id/edit" element={<EditToBareRedirect />} />
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
        {/* Phase B: `/edit` is gone — redirect old deep-links to the bare URL. */}
        <Route path="/compare/all/:id/edit" element={<EditToBareRedirect />} />
        {/* New permalink shapes — all render PageBody, which inspects Zustand
            to decide which page to mount. The URL-sync hooks dispatch state
            based on the matched route, so PageBody only needs to read the
            already-populated Zustand. */}
        <Route path="/" element={<PageBody />} />
        <Route path="/index" element={<PageBody />} />
        <Route path="/index/:experiment" element={<PageBody />} />
        <Route path="/index/:experiment/:sample" element={<PageBody />} />
        <Route path="/inspect" element={<PageBody />} />
        <Route path="/inspect/:experiment" element={<PageBody />} />
        <Route path="/inspect/:experiment/:sample" element={<PageBody />} />

        <Route path="*" element={<PageBody />} />  {/* stale fallback */}
      </Routes>
      <NavModal />
    </div>
  );
}
