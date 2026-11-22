import type { IngestProblems, ProblemExposure } from "../api";
import { Overlay } from "./Overlay";

interface Props {
  problems: IngestProblems | undefined;
  onPick: (sampleId: number) => void;
  onAddSetups: () => void;
  onClose: () => void;
}

function reason(p: ProblemExposure): string {
  return p.problem === "no_trace" ? "No trace yet. It is picked up when its _tot.dat appears." : p.reason;
}

export function ProblemList({ problems, onPick, onAddSetups, onClose }: Props): JSX.Element {
  const rows = problems?.exposures ?? [];
  const folders = problems?.setups ?? [];
  return (
    <Overlay label="Ingest problems" onClose={onClose} testId="workspace-problems">
      <div className="flex max-h-[70vh] flex-col gap-3 overflow-y-auto px-5 py-4 text-body">
        <h2 className="text-label text-ink-soft">Ingest problems</h2>
        {problems?.ingest_error != null && (
          <div data-testid="problems-ingest-error" className="flex flex-col gap-1">
            <span className="text-error">The last ingest failed.</span>
            <span className="whitespace-pre-wrap break-words text-caption text-ink-soft">{problems.ingest_error}</span>
          </div>
        )}
        {folders.length > 0 && (
          <div data-testid="problems-setups" className="flex flex-col gap-1">
            <span>{folders.length === 1 ? "A setup folder was not added" : `${folders.length} setup folders were not added`}: {folders.map((f) => f.folder).join(", ")}.</span>
            <span className="text-caption text-ink-soft">Their file layout was not recognised, so a rescan skips them. Add them with their file patterns typed in.</span>
            <button type="button" data-testid="problems-add-setups" onClick={onAddSetups} className="self-start font-semibold text-accent hover:underline">Add setup folders</button>
          </div>
        )}
        {rows.length === 0 && folders.length === 0 && problems?.ingest_error == null && <p className="text-ink-soft">No problems.</p>}
        {rows.length > 0 && (
          <ul className="flex flex-col">
            {rows.map((p) => (
              <li key={p.exposure_id} data-testid="problem-row" className="border-b border-hair py-2 last:border-b-0">
                <button
                  type="button"
                  disabled={p.sample_id === null}
                  onClick={() => { if (p.sample_id !== null) onPick(p.sample_id); }}
                  className="flex w-full flex-col items-start gap-0.5 text-left hover:text-accent"
                >
                  <span>{p.filename ?? "Unnamed frame"}{p.sample_name ? ` · ${p.sample_name}` : ""}</span>
                  <span className="text-caption text-ink-soft">{reason(p)}</span>
                </button>
              </li>
            ))}
          </ul>
        )}
      </div>
    </Overlay>
  );
}
