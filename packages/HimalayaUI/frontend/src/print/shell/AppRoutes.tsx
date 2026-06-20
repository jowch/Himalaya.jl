import { Routes, Route, Navigate, useParams, useNavigate } from "react-router-dom";
import { useAppState } from "../../state";
import { useGlobalShortcuts } from "../../hooks/useGlobalShortcuts";
import { CorpusShell } from "./CorpusShell";
import { ExperimentShell } from "./ExperimentShell";
import { SamplesPage } from "../pages/SamplesPage";
import { LoupePage } from "../pages/LoupePage";
import { FocusPage } from "../pages/FocusPage";
import { SeriesFolioPage } from "../pages/SeriesFolioPage";
import { SeriesBuilderPage } from "../pages/SeriesBuilderPage";
import { SeriesScopingPage } from "../pages/SeriesScopingPage";
import { ExperimentsHomePage } from "../pages/ExperimentsHomePage";
import { NewExperimentPage } from "../pages/NewExperimentPage";
import { ExperimentCorpusPage } from "../pages/ExperimentCorpusPage";
import { ConfigurationPage } from "../pages/ConfigurationPage";
import { GroupingReviewPage } from "../components/GroupingReviewPage";
import { IndexSlugRedirect } from "./IndexSlugRedirect";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";
import { useStateFromUrl } from "../../hooks/useStateFromUrl";

/** Thin wrapper that reads :id from the route and passes it to GroupingReviewPage. */
function GroupingReviewRoute(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const experimentId = id ? Number(id) : 0;
  return (
    <GroupingReviewPage
      experimentId={experimentId}
      onBack={() => navigate("../corpus")}
    />
  );
}

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
/**
 * SeriesBuilderRoute — keys the builder on its `:id` so each series mounts a
 * fresh instance.
 *
 * BU-MODESWITCH-LEAK: the builder holds id-scoped LOCAL state (the view
 * offset/scale, the draft undo/redo stacks, the Confirm-chain refs). React
 * reuses the same element instance across `/series/:id` param changes, so
 * without a per-id key those bleed from one series into the next (series 2's
 * view offset and undo history showing up on series 1). The global draft slot is
 * intentionally recoverable on return; this only resets the ephemeral per-visit
 * state that should never have crossed a series boundary.
 */
function SeriesBuilderRoute(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  return <SeriesBuilderPage key={id} />;
}

export function AppRoutes(): JSX.Element {
  // R0a (#221): the theme-class effect is gone. "The Print" is the single
  // identity defined statically in styles.css `@theme`; there is no
  // `theme-light` class to toggle on <html>.

  // Global keyboard shortcuts — hoisted above the shell so they work app-wide.
  // KEYS-LIB step 2: this now binds only the find/jump chord (`/`, `⌘K`); the
  // sample step `[`/`]` and every other gesture are surface-owned through the
  // shortcut library (print/shell/shortcuts.ts), so the per-experiment samples
  // list this used to thread for the old `,`/`.` stepper is no longer needed.
  useGlobalShortcuts();

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
        {/* I3.5a (#175): series builder — read-only visual surface. Keyed on
            :id (SeriesBuilderRoute) so per-series visits don't leak local state. */}
        <Route path="/series/:id" element={<SeriesBuilderRoute />} />
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
      <Route path="/" element={<Navigate to="/experiments" replace />} />
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
      {/* Ingestion redesign (spec §7/§9.6): the experiments tree sits OUTSIDE
          CorpusShell so ExperimentShell's own chrome never stacks on
          CorpusTopbar. */}
      <Route path="/experiments" element={<ExperimentsHomePage />} />
      <Route path="/experiments/new" element={<NewExperimentPage />} />
      <Route path="/experiments/:id" element={<ExperimentShell />}>
        <Route index element={<Navigate to="corpus" replace />} />
        <Route path="corpus" element={<ExperimentCorpusPage />} />
        <Route path="config" element={<ConfigurationPage />} />
        {/* E2: GroupingReviewPage mounts here (Task 20). */}
        <Route path="grouping" element={<GroupingReviewRoute />} />
      </Route>
    </Routes>
  );
}
