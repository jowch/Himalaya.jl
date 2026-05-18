import { useEffect } from "react";
import { Routes, Route, Navigate, useLocation } from "react-router-dom";
import { useAppState } from "../state";
import { useSamples } from "../queries";
import { useGlobalShortcuts } from "../hooks/useGlobalShortcuts";
import { CorpusShell } from "./CorpusShell";
import { SamplesPage } from "../pages/SamplesPage";
import { AppShell } from "./AppShell";
import { IndexPage } from "../pages/IndexPage";
import { InspectPage } from "../pages/InspectPage";
import { Compare } from "../pages/Compare";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";

/**
 * PageBody — the legacy Index/Inspect body switcher, driven by Zustand
 * `activePage`. Compare URLs are matched by their explicit <Route> entries
 * below, so activePage === "compare" never reaches here.
 */
function PageBody(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const resolving = useAppState((s) => s.resolving);
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  if (resolving) return <ResolvingFallback />;
  if (staleUrlContext !== null) {
    return <StaleUrlPage staleUrlContext={staleUrlContext} />;
  }
  if (activePage === "index") return <IndexPage />;
  if (activePage === "inspect") return <InspectPage />;
  return <></>;
}

/**
 * EditToBareRedirect — the `/edit` URL segment is gone (Compare UX Phase B).
 * Old `/edit` deep-links resolve by redirecting to the bare path. Reads
 * `useLocation()` (the router store), not `window.location` — under
 * MemoryRouter the latter stays at "/".
 */
function EditToBareRedirect(): JSX.Element {
  const loc = useLocation();
  const pathname = loc.pathname.replace(/\/edit\/?$/, "");
  return (
    <Navigate
      to={{ pathname, search: loc.search, hash: loc.hash }}
      replace
    />
  );
}

/**
 * AppRoutes — the single hoisted top-level <Routes> table, plus the shared
 * root effects (theme, global shortcuts) that sit above both shell bodies.
 *
 * Two pathless layout routes: <CorpusShell> for new corpus surfaces and
 * <AppShell> for the legacy Index/Inspect/Compare surfaces. Later redesign
 * issues register their route slot under the corpus layout route (#161 the
 * loupe, #179 the focus workspace).
 */
export function AppRoutes(): JSX.Element {
  const theme = useAppState((s) => s.theme);
  const experimentId = useAppState((s) => s.activeExperimentId);

  // Apply theme to <html> so CSS can key off `html.theme-light`. Hoisted
  // here (above both shell bodies) so the theme works under either shell.
  useEffect(() => {
    document.documentElement.className =
      theme === "light" ? "theme-light" : "";
    return () => {
      document.documentElement.className = ""; // defensive: StrictMode double-invoke / symmetry; not load-bearing
    };
  }, [theme]);

  // Global keyboard shortcuts — hoisted above both shell bodies so they work
  // under either. `useSamples(experimentId ?? 0)` matches the prior AppShell
  // call site; the `,`/`.` sample-step shortcut needs the active
  // experiment's samples. These shortcuts are now genuinely app-wide: they
  // fire under the corpus shell too (e.g. on /samples). #160 (contact sheet)
  // should be aware of this when landing — shortcuts may need guarding there.
  const samplesQ = useSamples(experimentId ?? 0);
  useGlobalShortcuts(experimentId === undefined ? undefined : samplesQ.data);

  return (
    <Routes>
      <Route element={<CorpusShell />}>
        <Route path="/samples" element={<SamplesPage />} />
      </Route>
      <Route element={<AppShell />}>
        <Route path="/" element={<PageBody />} />
        <Route path="/index" element={<PageBody />} />
        <Route path="/index/:experiment" element={<PageBody />} />
        <Route path="/index/:experiment/:sample" element={<PageBody />} />
        <Route path="/inspect" element={<PageBody />} />
        <Route path="/inspect/:experiment" element={<PageBody />} />
        <Route path="/inspect/:experiment/:sample" element={<PageBody />} />
        <Route path="/experiments/:eid/compare" element={<Compare />} />
        <Route path="/experiments/:eid/compare/new" element={<Compare />} />
        <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
        <Route
          path="/experiments/:eid/compare/:id/edit"
          element={<EditToBareRedirect />}
        />
        <Route path="/compare/all" element={<Compare />} />
        <Route path="/compare/all/new" element={<Compare />} />
        <Route path="/compare/all/:id" element={<Compare />} />
        <Route path="/compare/all/:id/edit" element={<EditToBareRedirect />} />
        <Route path="*" element={<PageBody />} />
      </Route>
    </Routes>
  );
}
