import { Overlay } from "./Overlay";

interface Props {
  name: string;
  samples: number;
  removing: boolean;
  error: string | null;
  onRemove: () => void;
  onClose: () => void;
}

export function RemoveExperimentDialog({ name, samples, removing, error, onRemove, onClose }: Props): JSX.Element {
  return (
    <Overlay label="Remove experiment" onClose={onClose} testId="workspace-remove-experiment">
      <div className="flex flex-col gap-3 px-5 py-4 text-body">
        <h2 className="text-label text-ink">Remove {name}?</h2>
        <p className="text-ink-soft">
          Its {samples.toLocaleString("en-US")} {samples === 1 ? "sample" : "samples"} and their phase calls, peaks, labels
          and bridges are deleted from Himalaya. The files in the data folder are not touched, so you can add the
          folder again from New experiment.
        </p>
        {error !== null && <p role="alert" className="text-error">{error}</p>}
        <div className="flex justify-end gap-2">
          <button type="button" onClick={onClose} className="rounded-sm border border-hair px-3 py-1 hover:bg-paper-sunk">
            Cancel
          </button>
          <button type="button" onClick={onRemove} disabled={removing}
            className="rounded-sm border border-error px-3 py-1 text-error hover:bg-paper-sunk disabled:opacity-50">
            {removing ? "Removing…" : "Remove experiment"}
          </button>
        </div>
      </div>
    </Overlay>
  );
}
