interface Props {
  progress: string | undefined;
  ingestFailed: boolean;
  waiting: number;
  failed: number;
  unrecognized: number;
  onDetails: () => void;
}

const plural = (n: number, one: string, many: string): string => `${n.toLocaleString("en-US")} ${n === 1 ? one : many}`;

export function ExperimentNotice({ progress, ingestFailed, waiting, failed, unrecognized, onDetails }: Props): JSX.Element | null {
  const parts = [
    ...(progress !== undefined ? [progress] : []),
    ...(waiting > 0 ? [`${plural(waiting, "frame", "frames")} waiting for a trace`] : []),
    ...(failed > 0 ? [`${plural(failed, "analysis", "analyses")} failed`] : []),
    ...(unrecognized > 0 ? [`${plural(unrecognized, "setup folder", "setup folders")} not recognised`] : []),
  ];
  if (parts.length === 0 && !ingestFailed) return null;
  const hasDetails = ingestFailed || waiting > 0 || failed > 0 || unrecognized > 0;
  return (
    <div data-testid="workspace-notice" role="status" className="flex items-baseline justify-between gap-2 border-b border-hair px-2.5 py-2 text-caption">
      <span className="flex flex-col gap-0.5">
        {ingestFailed && <span data-testid="workspace-ingest-failed" className="text-error">The last ingest failed.</span>}
        {parts.length > 0 && <span data-testid="workspace-ingest-status">{parts.join(" · ")}</span>}
      </span>
      {hasDetails && (
        <button type="button" onClick={onDetails} className="whitespace-nowrap font-semibold text-accent hover:underline">
          Details
        </button>
      )}
    </div>
  );
}
