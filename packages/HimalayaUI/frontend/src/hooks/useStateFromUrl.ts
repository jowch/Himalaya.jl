import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { parseLocation } from "../lib/url/parseLocation";

// Spec §4.2 — URL → Zustand. Reads `useLocation()` so popstate AND
// useNavigate both flow through the hook.
//
// I4.4 (#181) + I1.7 (#163) + I3.6 (#177): the Index, Inspect, and Compare
// surfaces are all retired. `parseLocation` no longer returns a slug-bearing
// `index`/`inspect`/`compare` kind, so this hook no longer resolves slug URLs
// into Zustand active ids — those URLs are redirected at the router (Index →
// /samples or /sample/:id; Inspect → /samples; Compare → /series). What
// remains:
//   - `root` (/) → redirect to the corpus contact sheet (/samples) per §4.1.
//   - `stale`    → park the user on StaleUrlPage.

export function useStateFromUrl(): void {
  const location = useLocation();
  const navigate = useNavigate();

  useEffect(() => {
    const parsed = parseLocation(location.pathname, location.search);

    if (parsed.kind === "stale") {
      useAppState.getState().setStaleUnknownPath(parsed.raw);
      return;
    }

    // parsed.kind === "root". §4.1: the app home is the corpus contact sheet.
    // Index is retired (#181) and Compare is URL-owned, so a cold `/` has
    // nothing to reflect from Zustand — send the user to the corpus.
    //
    // NOTE: this arm is effectively test-only in production. Bare `/` is
    // intercepted by the standalone `<Route path="/">` redirect in AppRoutes
    // (registered OUTSIDE the AppShell layout that hosts this hook), so a real
    // browser never mounts AppShell at `/` and `parseLocation` never returns
    // `root` here. It's exercised by useStateFromUrl's unit tests, and kept as
    // a correct belt-and-suspenders fallback should the route ever change.
    navigate("/samples", { replace: true });
  }, [location.pathname, location.search, navigate]);
}
