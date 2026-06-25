import { Routes, Route, Navigate, useParams, useNavigate } from "react-router-dom";
import { Outlet } from "react-router-dom";
import { useAppState } from "../../state";
import { useGlobalShortcuts } from "../../hooks/useGlobalShortcuts";
import { TopNav } from "./TopNav";
import { ExperimentShell } from "./ExperimentShell";
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
import { InteractionDock } from "../interaction/InteractionDock";
import { useKeyboardLayer } from "../interaction/useKeyboardLayer";

/** Thin wrapper that reads :id from the route and passes it to GroupingReviewPage. */
function GroupingReviewRoute(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const experimentId = id ? Number(id) : 0;
  return (
    <GroupingReviewPage
      // Remount per experiment so ephemeral selection/cursor state can't carry a
      // stale multi-selection across an experiment switch (mirrors SeriesBuilderRoute).
      key={experimentId}
      experimentId={experimentId}
      // Absolute target: `/grouping` is a top-level takeover route (not nested
      // under the experiment shell), so a relative `../corpus` resolves to bare
      // `/corpus` (a 404). Confirm-groups and Back both land on this
      // experiment's corpus.
      onBack={() => navigate(`/experiments/${experimentId}/corpus`)}
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

/**
 * AppRoutes — the single hoisted top-level <Routes> table, plus the shared
 * root effects (global shortcuts).
 *
 * T3.2 (App Shell Unification): one layout route, <AppShell>, hosts every
 * surface — the corpus pages, the experiment tree, Focus/Loupe, Series, AND
 * the `*` stale catch-all (PageBody). The legacy shell trio was
 * deleted in T3.2. Pre-shell redirects (`/`, `/index*`,
 * `/compare*`, `/samples`) sit outside the layout route so no chrome mounts
 * under them and races the redirect.
 *
 * Focus/Loupe stay flat at `/sample/:id` and `/sample/:id/loupe` — NO resolver,
 * NO experiment-scoped nesting (a resolver would always flash a loading screen).
 * The chrome reads the experiment from the loaded sample's experiment_id.
 */
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
      {/* T3.2: single AppShell wraps every surface. */}
      <Route element={<AppShell />}>
        {/* T3.2: flat loupe route — /sample/:sampleId/loupe (matches gotoSample + loupeHref). */}
        <Route path="/sample/:sampleId/loupe" element={<LoupePage />} />
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
        {/* I1.7 (#163): Inspect retired. Old /inspect* deep-links land on
            /experiments (the new home). Splat covers /inspect, /inspect/:exp,
            /inspect/:exp/:sample. */}
        <Route path="/inspect/*" element={<Navigate to="/experiments" replace />} />
        {/* Ingestion redesign: experiments tree (E1). ExperimentShell owns its
            own header/tabs per spec §3.2 (page content, not a separate chrome). */}
        <Route path="/experiments" element={<ExperimentsHomePage />} />
        <Route path="/experiments/new" element={<NewExperimentPage />} />
        {/* T4.0: draft Configuration — first-run mode, no :id yet. */}
        <Route path="/experiments/new/config" element={<ConfigurationPage />} />
        {/* The combined scan + grouping-review surface (p1-grouping) is a
            TAKEOVER: it owns its full header, so it mounts OUTSIDE ExperimentShell
            (no redundant experiment header). Declared before the layout route. */}
        <Route path="/experiments/:id/grouping" element={<GroupingReviewRoute />} />
        <Route path="/experiments/:id" element={<ExperimentShell />}>
          <Route index element={<Navigate to="corpus" replace />} />
          <Route path="corpus" element={<ExperimentCorpusPage />} />
          <Route path="config" element={<ConfigurationPage />} />
        </Route>
        {/* The `*` stale catch-all mounts PageBody which runs the URL→stale
            classifier and renders StaleUrlPage or ResolvingFallback. */}
        <Route path="*" element={<PageBody />} />
      </Route>
      {/* Pre-shell redirects — sit OUTSIDE the layout route so no chrome
          mounts under them and races the redirect. */}
      {/* I4.4 (#181): Index retired. Bare `/` and the sampleless legacy Index
          URLs land on the experiments home; a slug-bearing `/index/:exp/:sample`
          resolves to the focus workspace via IndexSlugRedirect. */}
      <Route path="/" element={<Navigate to="/experiments" replace />} />
      <Route path="/index" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment" element={<Navigate to="/experiments" replace />} />
      <Route path="/index/:experiment/:sample" element={<IndexSlugRedirect />} />
      {/* T3.2: /samples redirects to /experiments (SamplesPage retired from /samples). */}
      <Route path="/samples" element={<Navigate to="/experiments" replace />} />
      {/* I3.6 (#177): Compare retired. The series stage replaces it; all
          `/compare*` deep-links redirect to the series folio. */}
      <Route path="/experiments/:eid/compare/*" element={<Navigate to="/series" replace />} />
      <Route path="/compare/all/*" element={<Navigate to="/series" replace />} />
    </Routes>
  );
}
