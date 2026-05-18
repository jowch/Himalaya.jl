import { useEffect } from "react";
import { Outlet, useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { useStateFromUrl } from "../hooks/useStateFromUrl";
import { useUrlFromState } from "../hooks/useUrlFromState";

/**
 * AppShell — the legacy layout: the grain background, the app header, the
 * page-nav rocker, and an <Outlet/> for the matched legacy route (Index /
 * Inspect / Compare). Registered as a pathless layout route in `AppRoutes`.
 *
 * Phase 5 (#182) deletes this; until then it is the migration rollback path.
 *
 * Routing model: AppShell is the element of the *legacy* layout route, so it
 * mounts only when a legacy URL matches. The URL↔Zustand sync hooks live
 * here (not at the app root) precisely so they never run on a corpus route —
 * the new shell owns its own URL and the legacy stale-path logic cannot
 * strand the user there.
 */
export function AppShell(): JSX.Element {
  // URL ↔ Zustand sync, relocated from App.tsx. Order matters:
  // useStateFromUrl populates Zustand before useUrlFromState reflects it
  // back — the `resolving`-flag handshake between the two depends on
  // useStateFromUrl's effect running first.
  useStateFromUrl();
  useUrlFromState();

  const experimentId = useAppState((s) => s.activeExperimentId);
  const setActivePage = useAppState((s) => s.setActivePage);
  const activePage = useAppState((s) => s.activePage);
  const location = useLocation();
  const navigate = useNavigate();

  // Sync URL → Zustand activePage. When the URL is a compare path, mark the
  // page tab as "compare". On other paths we leave activePage alone.
  const onComparePath =
    location.pathname.startsWith("/compare") ||
    /^\/experiments\/\d+\/compare(\/|$)/.test(location.pathname);
  useEffect(() => {
    if (onComparePath) setActivePage("compare");
  }, [onComparePath, setActivePage]);

  // Symmetric: when activePage is "compare" but the URL isn't on a compare
  // path, navigate so the URL-routed Compare page mounts. Without this, a
  // reload at "/" with persisted activePage='compare' renders the rocker but
  // no page body (issue #77).
  useEffect(() => {
    if (activePage !== "compare") return;
    if (onComparePath) return;
    const url =
      experimentId !== undefined
        ? `/experiments/${experimentId}/compare`
        : "/compare/all";
    navigate(url, { replace: true });
  }, [activePage, onComparePath, experimentId, navigate]);

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full max-w-[1600px] mx-auto flex flex-col min-h-0 relative"
      // --chrome-h: AppHeader (h-11 = 44px) + TabRocker row (~40px).
      style={{ "--chrome-h": "84px" } as React.CSSProperties}
    >
      <AppHeader />
      <div className="shrink-0 flex justify-center pt-1 pb-2">
        <TabRocker
          experimentId={experimentId}
          onNavigateAway={(target) => {
            // Leaving Compare → return to "/" so the legacy body renders the
            // chosen page.
            if (onComparePath && target !== "compare") navigate("/");
          }}
        />
      </div>
      <Outlet />
      <NavModal />
    </div>
  );
}
