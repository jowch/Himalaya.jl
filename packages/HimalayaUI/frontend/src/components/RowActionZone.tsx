interface Props {
  onOverflow: () => void;
}

export function RowActionZone({ onOverflow }: Props): JSX.Element {
  return (
    <div className="inline-flex items-center gap-1 text-ink-faint">
      <button
        type="button"
        data-testid="row-action-overflow"
        data-interactable="button"
        onClick={(e) => { e.stopPropagation(); onOverflow(); }}
        className="px-1 hover:text-ink hover:bg-paper-sunk rounded"
        title="Row actions"
      >
        ⋯
      </button>
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
