import { useEffect } from "react";
import { Routes, Route, Navigate, useLocation } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
import { useAppState } from "../state";
import { queryKeys } from "../queries";
import * as api from "../api";
import { useGlobalShortcuts } from "../hooks/useGlobalShortcuts";
import { CorpusShell } from "./CorpusShell";
import { SamplesPage } from "../pages/SamplesPage";
import { LoupePage } from "../pages/LoupePage";
import { FocusWorkspacePage } from "../pages/FocusWorkspacePage";
import { SeriesFolioPage } from "../pages/SeriesFolioPage";
import { SeriesBuilderPage } from "../pages/SeriesBuilderPage";
import { SeriesScopingPage } from "../pages/SeriesScopingPage";
import { AppShell } from "./AppShell";
import { IndexPage } from "../pages/IndexPage";
import { Compare } from "../pages/Compare";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";

/**
 * PageBody — the legacy Index body, driven by Zustand `activePage`. Compare
 * URLs are matched by their explicit <Route> entries below, so
 * activePage === "compare" never reaches here. (Inspect retired in #163.)
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
 * <AppShell> for the legacy Index/Compare surfaces. Later redesign
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
  // under either; the `,`/`.` sample-step shortcut needs the active
  // experiment's samples. These shortcuts are now genuinely app-wide: they
  // fire under the corpus shell too (e.g. on /samples). #160 (contact sheet)
  // should be aware of this when landing — shortcuts may need guarding there.
  //
  // Gated on an active experiment via `useQuery` directly (the precedent in
  // useUrlFromState.ts) — `useSamples` exposes no `enabled` option, and
  // without the gate this fires GET /api/samples?experiment=0 on every cold
  // mount before an experiment is picked.
  const samplesQuery = useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId as number),
    enabled: experimentId !== undefined,
  });
  useGlobalShortcuts(
    experimentId === undefined ? undefined : samplesQuery.data,
  );

  return (
    <Routes>
      <Route element={<CorpusShell />}>
        <Route path="/samples" element={<SamplesPage />} />
        <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        {/* I4.1 (#178): focus workspace. I4.4 (#181) redirects /index* here. */}
        <Route path="/sample/:sampleId" element={<FocusWorkspacePage />} />
        {/* I3.3 (#173): series folio — corpus-wide masonry of saved series. */}
        <Route path="/series" element={<SeriesFolioPage />} />
        {/* I3.4 (#174): series scoping — the confirm-and-build gate that writes
            scoping sample_tags, then lands on the folio. Static /series/new
            ranks above the dynamic /series/:id (react-router v6 specificity). */}
        <Route path="/series/new" element={<SeriesScopingPage />} />
        {/* I3.5a (#175): series builder — read-only visual surface. */}
        <Route path="/series/:id" element={<SeriesBuilderPage />} />
        {/* I1.7 (#163): Inspect retired. Old /inspect* deep-links land on the
            contact sheet. Splat covers /inspect, /inspect/:exp, /inspect/:exp/:sample. */}
        <Route path="/inspect/*" element={<Navigate to="/samples" replace />} />
      </Route>
      <Route element={<AppShell />}>
        <Route path="/" element={<PageBody />} />
        <Route path="/index" element={<PageBody />} />
        <Route path="/index/:experiment" element={<PageBody />} />
        <Route path="/index/:experiment/:sample" element={<PageBody />} />
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
