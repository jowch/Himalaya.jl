interface Props {
  rescan: string | undefined;
  failed: boolean;
  flagged: number;
  onReview: () => void;
  onFailure: () => void;
}

export function ExperimentNotice({ rescan, failed, flagged, onReview, onFailure }: Props): JSX.Element | null {
  if (rescan === undefined && !failed && flagged === 0) return null;
  return (
    <div data-testid="workspace-notice" role="status" className="flex flex-col gap-1 border-b border-hair px-2.5 py-2 text-caption">
      {rescan !== undefined && <span data-testid="workspace-rescan">{rescan}</span>}
      {failed && (
        <span data-testid="workspace-rescan-failed" className="flex items-baseline justify-between gap-2">
          <span className="text-error">The last rescan failed.</span>
          <button type="button" onClick={onFailure} className="whitespace-nowrap font-semibold text-accent hover:underline">
            Details
          </button>
        </span>
      )}
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
