/**
 * RowActionZone — the drag cue (⋮⋮) on a compare/series member row. M2 removed
 * the sibling `⋯` overflow button: it only toggled a `data-overflow-open`
 * attribute that nothing consumed, opening no menu (a no-op affordance). The
 * remaining `⋮⋮` is inert signage for the real reorder grip; the caller renders
 * it only in edit mode, where that grip actually exists.
 */
export function RowActionZone(): JSX.Element {
  return (
    <div className="inline-flex items-center text-ink-faint">
      <span
        data-testid="row-action-drag-cue"
        aria-hidden="true"
        className="px-1 select-none text-ink-faint cursor-grab"
      >
        ⋮⋮
      </span>
    </div>
  );
}
