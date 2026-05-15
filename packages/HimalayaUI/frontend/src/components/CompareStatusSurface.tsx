import { useEffect, useState } from "react";

interface NeedsReview { onReanalyze: () => void; }
// `previousHash` is the hash the user last acknowledged; the banner
// surfaces while it differs from the current `content_hash`.
// `onAcknowledge` rebases the acknowledged hash to the current value
// (see `useStaleAgainstHash` in C-15 Step 4b).
interface ServerUpdate { previousHash: string; onAcknowledge: () => void; }

interface Props {
  needsReview: NeedsReview | null;
  serverUpdate: ServerUpdate | null;
  savedAt: number | null;
}

export function CompareStatusSurface(p: Props): JSX.Element | null {
  // `savedAt` is a wall-clock timestamp; React does not re-render on the
  // passage of time, so the "Saved" pill needs an explicit timeout to
  // retire itself 2s after the save. Without the effect the pill would
  // stick forever once shown.
  // Initialise from the prop so a save that already happened shows the pill
  // on the very first render — avoids a one-frame flash of nothing before
  // the effect below fires.
  const [showSaved, setShowSaved] = useState<boolean>(
    () => p.savedAt !== null && Date.now() - p.savedAt < 2000,
  );
  useEffect(() => {
    if (p.savedAt === null) {
      setShowSaved(false);
      return;
    }
    const elapsed = Date.now() - p.savedAt;
    if (elapsed >= 2000) {
      setShowSaved(false);
      return;
    }
    setShowSaved(true);
    const handle = setTimeout(() => setShowSaved(false), 2000 - elapsed);
    return () => clearTimeout(handle);
  }, [p.savedAt]);

  if (!p.needsReview && !p.serverUpdate && !showSaved) return null;
  return (
    <div data-testid="compare-status-surface" className="flex flex-col gap-1">
      {p.needsReview && (
        <div className="px-3 py-2 rounded border border-warning bg-warning/10 text-sm text-fg flex items-center gap-2">
          <span aria-hidden="true">⚠</span>
          <span>
            Needs review — analysis changed since last submit.
          </span>
          <button
            type="button"
            data-testid="compare-status-resnapshot"
            data-interactable="button"
            onClick={p.needsReview.onReanalyze}
            className="ml-auto px-2 py-0.5 rounded border border-border hover:bg-bg-hover text-xs"
          >
            Re-snapshot
          </button>
        </div>
      )}
      {p.serverUpdate && (
        <div
          data-testid="compare-status-server-update"
          className="px-3 py-2 rounded border border-border bg-bg-subtle text-sm text-fg flex items-center gap-2"
        >
          <span>Server-side updated since you last viewed — save may conflict.</span>
          <button
            type="button"
            data-testid="compare-status-acknowledge"
            data-interactable="button"
            onClick={p.serverUpdate.onAcknowledge}
            className="ml-auto px-2 py-0.5 rounded border border-border hover:bg-bg-hover text-xs"
          >
            Acknowledge
          </button>
        </div>
      )}
      {showSaved && (
        <div className="px-3 py-2 rounded border border-success bg-success/10 text-sm text-fg">
          Saved
        </div>
      )}
    </div>
  );
}
