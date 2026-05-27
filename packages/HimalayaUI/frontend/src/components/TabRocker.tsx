import type { PageId } from "../state";

// I3.6 (#177): Compare retired (Index #181, Inspect #163 before it) — there are
// no legacy page tabs left, so the rocker renders nothing. The component (and
// the whole dual-nav model) is deleted in I5.1; kept as an inert stub until
// then so AppShell's layout doesn't need restructuring twice.
const TABS: readonly { id: PageId; label: string }[] = [];

interface Props {
  /** Active experiment id — unused now that no tab builds a legacy URL.
   *  Retained until I5.1 deletes the shell + rocker together. */
  experimentId?: number | undefined;
  /** Called when the user leaves the (now-empty) rocker. Unused while TABS is
   *  empty; retained for the I5.1 deletion seam. */
  onNavigateAway?: (target: PageId) => void;
}

/**
 * TabRocker — pill-style segmented control for app-level page switching.
 * With every legacy surface retired (#163/#181/#177) there are no tabs to
 * render; this is an inert stub that I5.1 deletes along with the dual-nav
 * model. Props are kept so AppShell needn't change twice.
 */
export function TabRocker({ experimentId, onNavigateAway }: Props): JSX.Element {
  void experimentId;
  void onNavigateAway;
  return (
    <div
      role="tablist"
      data-testid="tab-rocker"
      className="inline-flex items-center gap-0.5 p-0.5
                 bg-bg-elevated border border-border rounded-full"
    >
      {TABS.map((t) => (
        <button
          key={t.id}
          role="tab"
          data-testid={`tab-${t.id}`}
          className="px-3.5 py-1 rounded-full font-sans text-sm font-medium text-fg-muted"
        >
          {t.label}
        </button>
      ))}
    </div>
  );
}
