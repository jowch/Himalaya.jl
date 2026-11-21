interface Props {
  rescan: string | undefined;
  flagged: number;
  onReview: () => void;
}

export function ExperimentNotice({ rescan, flagged, onReview }: Props): JSX.Element | null {
  if (rescan === undefined && flagged === 0) return null;
  return (
    <div data-testid="workspace-notice" role="status" className="flex flex-col gap-1 border-b border-hair px-2.5 py-2 text-caption">
      {rescan !== undefined && <span data-testid="workspace-rescan">{rescan}</span>}
      {flagged > 0 && (
        <span className="flex items-baseline justify-between gap-2">
          <span>{flagged === 1 ? "1 sample needs" : `${flagged} samples need`} a grouping check.</span>
          <button type="button" onClick={onReview} className="whitespace-nowrap font-semibold text-accent hover:underline">
            Review
          </button>
        </span>
      )}
    </div>
  );
}
