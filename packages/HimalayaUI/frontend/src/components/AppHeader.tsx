import { UtilityCluster } from "./UtilityCluster";

/**
 * AppHeader — the global, app-scoped header row.
 *
 * Layout:
 *   [               empty               ] [ UtilityCluster (right) ]
 *
 * The page-nav TabRocker sits in its own row below this one (inside
 * AppShell) so the page title can claim the plot card's top strip
 * without being crowded by the page-nav control.
 */
export function AppHeader(): JSX.Element {
  return (
    <header
      data-testid="app-header"
      className="h-11 flex items-center justify-end px-3"
    >
      <UtilityCluster />
    </header>
  );
}
