import { Outlet, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { useStateFromUrl } from "../hooks/useStateFromUrl";
import { useUrlFromState } from "../hooks/useUrlFromState";

/**
 * AppShell — the legacy layout: the grain background, the app header, the
 * page-nav rocker, and an <Outlet/> for the surviving catch-all route. All
 * legacy surfaces are retired — Inspect (#163), Index (#181), Compare (#177) —
 * so AppShell now mounts only for the `*` catch-all (the stale/resolving body).
 *
 * Phase 5 (#182) deletes this whole shell + the `activePage` model; until then
 * it is the migration rollback path.
 *
 * Routing model: AppShell is the element of the *legacy* layout route, so it
 * mounts only when no other (corpus/redirect) route matches. The URL↔Zustand
 * sync hooks live here (not at the app root) precisely so they never run on a
 * corpus route — the new shell owns its own URL and the legacy stale-path
 * logic cannot strand the user there.
 */
export function AppShell(): JSX.Element {
  // URL ↔ Zustand sync, relocated from App.tsx. Order matters:
  // useStateFromUrl populates Zustand before useUrlFromState reflects it
  // back — the `resolving`-flag handshake between the two depends on
  // useStateFromUrl's effect running first.
  useStateFromUrl();
  useUrlFromState();

  const experimentId = useAppState((s) => s.activeExperimentId);
  const navigate = useNavigate();

  // I3.6 (#177): Compare is retired. The `/compare*` URLs redirect to /series
  // at the router (outside this shell), so AppShell never mounts on a compare
  // path and `activePage` is the inert "none" sentinel. The two former
  // URL↔activePage nav-bridge effects (and the `onComparePath` detection) are
  // gone with the Compare page; the dual-nav model retires entirely in I5.1.

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
          onNavigateAway={() => {
            // No legacy surface remains to navigate away from; return to the
            // corpus home. (TabRocker + this shell retire in I5.1.)
            navigate("/samples");
          }}
        />
      </div>
      <Outlet />
      <NavModal />
    </div>
  );
}
