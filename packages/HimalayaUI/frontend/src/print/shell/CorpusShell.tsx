import { Outlet } from "react-router-dom";
import { CorpusTopbar } from "./CorpusTopbar";

/**
 * CorpusShell — the layout-route element for the redesigned corpus-scoped
 * surfaces. Renders the topbar and an <Outlet/> for the matched child route.
 *
 * Registered as a pathless layout route in `AppRoutes`; later issues add
 * child routes under it (#161 the loupe, #179 the focus workspace).
 */
export function CorpusShell(): JSX.Element {
  return (
    <div
      data-testid="corpus-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink"
    >
      <CorpusTopbar />
      <main className="flex-1 min-h-0 overflow-auto">
        <Outlet />
      </main>
    </div>
  );
}
