import { useSyncActiveSampleFromRoute } from "../hooks/useSyncActiveSampleFromRoute";
import { FocusWorkspaceLayout } from "../components/FocusWorkspaceLayout";

/**
 * FocusWorkspacePage — the focus workspace at `/sample/:sampleId`. Mounted
 * under the CorpusShell layout route.
 *
 * The route param seeds Zustand `activeSampleId` via
 * `useSyncActiveSampleFromRoute` (I4.1 / #178). The body is the focus-workspace
 * layout (I4.2 / #179): trace hero, co-resident detector, phase rail, Notes
 * margin — all composed from existing index components reused unchanged.
 *
 * The q-link cross-highlight is I4.3 (#180); the `/index` cutover is I4.4 (#181).
 */
export function FocusWorkspacePage(): JSX.Element {
  useSyncActiveSampleFromRoute();
  return (
    <div data-testid="focus-workspace-page" className="flex min-h-0 flex-1 flex-col">
      <FocusWorkspaceLayout />
    </div>
  );
}
