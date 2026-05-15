interface Props {
  onOverflow: () => void;
}

export function RowActionZone({ onOverflow }: Props): JSX.Element {
  return (
    <div className="inline-flex items-center gap-1 text-fg-dim">
      <button
        type="button"
        data-testid="row-action-overflow"
        data-interactable="button"
        onClick={(e) => { e.stopPropagation(); onOverflow(); }}
        className="px-1 hover:text-fg hover:bg-bg-hover rounded"
        title="Row actions"
      >
        ⋯
      </button>
      <span
        data-testid="row-action-drag-cue"
        aria-hidden="true"
        className="px-1 select-none text-fg-dim cursor-grab"
      >
        ⋮⋮
      </span>
    </div>
  );
}
