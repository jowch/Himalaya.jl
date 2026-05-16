import { useLocation } from "react-router-dom";
import { useAppState } from "../state";
import { ComparePage } from "./ComparePage";
import { ComparePageEdit } from "./ComparePageEdit";

/**
 * Compare.tsx — unified single-mode shell for all compare routes.
 *
 * Compare UX C-11: this is a pure passthrough that delegates to
 * `ComparePageEdit` when a new/edit flow is active and `ComparePage` for
 * review. Subsequent Compare UX tasks (C-12 onward) progressively replace
 * the delegated bodies with the new title strip / toolbar / save pill, and
 * ultimately delete ComparePage.tsx + ComparePageEdit.tsx.
 *
 * The wrapper div carries `data-testid="compare-page"` (the canonical
 * testid for the compare surface) and `className="contents"` so it is
 * transparent to the flex/grid layout of its parent — it does not consume
 * a grid track in WorkspaceGrid's 3-column subgrid.
 */
export function Compare(): JSX.Element {
  const location = useLocation();
  const isNew = location.pathname.endsWith("/new");
  const hasDraft = useAppState((s) => s.activeDraft !== null);
  const Body = (isNew || hasDraft) ? ComparePageEdit : ComparePage;
  return (
    <div data-testid="compare-page" className="contents">
      <Body />
    </div>
  );
}
