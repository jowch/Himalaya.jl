import { useAppState, type PageId } from "../state";

const TABS: readonly { id: PageId; label: string }[] = [
  { id: "inspect", label: "Inspect" },
  { id: "index",   label: "Index"   },
  { id: "compare", label: "Compare" },
];

/**
 * TabRocker — pill-style segmented control for app-level page switching.
 * Shared across all pages; lives in AppHeader.
 */
export function TabRocker(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const setPage    = useAppState((s) => s.setActivePage);

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
            onClick={() => setPage(t.id)}
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
