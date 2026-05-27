import { useSyncActiveSampleFromRoute } from "../hooks/useSyncActiveSampleFromRoute";

/**
 * FocusWorkspacePage — the focus workspace at `/sample/:sampleId`. Mounted
 * under the CorpusShell layout route.
 *
 * I4.1 (#178): a wiring shim only. It seeds the Zustand `activeSampleId` from
 * the route param (via `useSyncActiveSampleFromRoute`) so the reused index
 * components work under the URL-routed surface, and renders a placeholder
 * body. I4.2 (#179) replaces the placeholder with the real trace-hero /
 * co-resident detector / phase-rail / Notes-drawer layout.
 */
export function FocusWorkspacePage(): JSX.Element {
  useSyncActiveSampleFromRoute();
  return (
    <div data-testid="focus-workspace-page" className="p-8 text-sm text-ink-faint">
      Focus workspace — layout lands in #179.
    </div>
  );
}
