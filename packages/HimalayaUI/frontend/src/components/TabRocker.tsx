import { useNavigate } from "react-router-dom";
import { useAppState, type PageId } from "../state";

// I4.4 (#181): Index retired (Inspect retired #163) — only Compare remains.
// The rocker (and the whole dual-nav model) is deleted in I5.1.
const TABS: readonly { id: PageId; label: string }[] = [
  { id: "compare", label: "Compare" },
];

interface Props {
  /**
   * Active experiment id (when set), used to build the Compare URL.
   * `/experiments/:eid/compare` if set, `/compare/all` otherwise.
   */
  experimentId?: number | undefined;
  /**
   * Called when the user clicks a tab other than the current Compare tab
   * while on a Compare URL — lets the AppShell navigate back to "/".
   */
  onNavigateAway?: (target: PageId) => void;
}

/**
 * TabRocker — pill-style segmented control for app-level page switching.
 * Shared across all pages; lives in AppHeader.
 *
 * Hybrid navigation: clicking "Compare" navigates to a URL (so :eid/:id
 * params survive reloads). The Index/Inspect tabs (which updated `activePage`
 * and let the parent shell route back to "/") are retired (#181 / #163).
 */
export function TabRocker({ experimentId, onNavigateAway }: Props): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const setPage    = useAppState((s) => s.setActivePage);
  const navigate   = useNavigate();

  const handleClick = (target: PageId): void => {
    setPage(target);
    if (target === "compare") {
      const url = experimentId !== undefined
        ? `/experiments/${experimentId}/compare`
        : "/compare/all";
      navigate(url);
    } else {
      onNavigateAway?.(target);
    }
  };

  return (
    <div
      role="tablist"
      data-testid="tab-rocker"
      className="inline-flex items-center gap-0.5 p-0.5
                 bg-bg-elevated border border-border rounded-full"
    >
      {TABS.map((t) => {
        const active = t.id === activePage;
        return (
          <button
            key={t.id}
            role="tab"
            aria-selected={active}
            data-testid={`tab-${t.id}`}
            data-active={active || undefined}
            onClick={() => handleClick(t.id)}
            className={
              "px-3.5 py-1 rounded-full font-sans text-sm font-medium " +
              (active
                ? "bg-accent/15 text-accent"
                : "text-fg-muted hover:text-fg")
            }
          >
            {t.label}
          </button>
        );
      })}
    </div>
  );
}
