import { Routes, Route, Navigate } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
import { useAppState } from "../../state";
import { queryKeys } from "../../queries";
import * as api from "../../api";
import { useGlobalShortcuts } from "../../hooks/useGlobalShortcuts";
import { CorpusShell } from "./CorpusShell";
import { SamplesPage } from "../pages/SamplesPage";
import { LoupePage } from "../pages/LoupePage";
import { FocusPage } from "../pages/FocusPage";
import { SeriesFolioPage } from "../pages/SeriesFolioPage";
import { SeriesBuilderPage } from "../pages/SeriesBuilderPage";
import { SeriesScopingPage } from "../pages/SeriesScopingPage";
import { IndexSlugRedirect } from "./IndexSlugRedirect";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";
import { useStateFromUrl } from "../../hooks/useStateFromUrl";

/**
 * PageBody — the `*` catch-all body. After I5.1 (#182) it mounts under the
 * single CorpusShell (the legacy AppShell + dual-nav model are retired). Index
 * (#181), Inspect (#163), and Compare (#177) are all retired, so it renders
 * the in-flight (`resolving`) and stale (`staleUrlContext`) URL states only.
 *
 * The URL→stale classifier (`useStateFromUrl`) is mounted HERE, inside the
 * catch-all element, NOT in any layout shell. PageBody renders only when no
 * other route matched, so the classifier runs only on a genuinely-unmatched
 * path — mounting it in CorpusShell would wrongly park a valid corpus route
 * (e.g. `/samples`, which `parseLocation` has no arm for) on StaleUrlPage.
 * This preserves the isolation the AppShell layout used to provide.
 *
 * Any path that reaches `*` is, by definition, one no other route matched —
 * `useStateFromUrl` classifies it as a stale unknown_path and sets
 * `staleUrlContext` inside its effect. On the first render that effect hasn't
 * run yet (so `staleUrlContext` is still null); we render the neutral
 * ResolvingFallback for that one frame rather than `<Navigate to="/samples">`.
 * Navigating away on the bare fallback raced the stale dispatch and bounced
 * typo'd/dead URLs off the "Page not found" view (the #181 regression).
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
 * AppRoutes — the single hoisted top-level <Routes> table, plus the shared
 * root effects (global shortcuts).
 *
 * I5.1 (#182): one layout route, <CorpusShell>, hosts every surface — the
 * corpus pages AND the `*` stale catch-all (PageBody). The legacy AppShell +
 * the dual-nav `activePage` model are retired. Pre-shell redirects (`/`,
 * `/index*`, `/compare*`) sit outside the layout route so no chrome mounts
 * under them and races the redirect.
 */
export function AppRoutes(): JSX.Element {
  const experimentId = useAppState((s) => s.activeExperimentId);

  // R0a (#221): the theme-class effect is gone. "The Print" is the single
  // identity defined statically in styles.css `@theme`; there is no
  // `theme-light` class to toggle on <html>.

  // Global keyboard shortcuts — hoisted above the shell so they work app-wide;
  // the `,`/`.` sample-step shortcut needs the active experiment's samples.
  // They fire under the corpus shell (e.g. on /samples). #160 (contact sheet)
  // should be aware of this when landing — shortcuts may need guarding there.
  //
  // Gated on an active experiment via `useQuery` directly — `useSamples`
  // exposes no `enabled` option, and without the gate this fires GET
  // /api/samples?experiment=0 on every cold mount before an experiment is picked.
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
        <Route path="/sample/:sampleId" element={<FocusPage />} />
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
        {/* I5.1 (#182): the `*` stale catch-all now lives under the single
            CorpusShell. PageBody mounts the URL→stale classifier; StaleUrlPage
            renders under CorpusTopbar. (Was the legacy AppShell layout route.) */}
        <Route path="*" element={<PageBody />} />
      </Route>
      {/* I4.4 (#181): Index retired. These redirects sit OUTSIDE the layout
          shell so no chrome mounts under them and races the redirect. Bare `/`
          and the sampleless legacy Index URLs land on the corpus contact sheet;
          a slug-bearing `/index/:exp/:sample` resolves to the focus workspace
          via IndexSlugRedirect (preserving old permalink deep-links). */}
      <Route path="/" element={<Navigate to="/samples" replace />} />
      <Route path="/index" element={<Navigate to="/samples" replace />} />
      <Route path="/index/:experiment" element={<Navigate to="/samples" replace />} />
      <Route path="/index/:experiment/:sample" element={<IndexSlugRedirect />} />
      {/* I3.6 (#177): Compare retired. The series stage replaces it; all
          `/compare*` deep-links (both the experiment-scoped and the global
          `/compare/all` roots, incl. `/new`, `/:id`, `/:id/edit`) redirect to
          the series folio. Placed OUTSIDE the layout shell for the same reason
          as the Index redirects above. The `comparison*` tables + dispatcher
          branches stay forever for event replay; only the UI + routes are gone. */}
      <Route path="/experiments/:eid/compare/*" element={<Navigate to="/series" replace />} />
      <Route path="/compare/all/*" element={<Navigate to="/series" replace />} />
    </Routes>
  );
}
