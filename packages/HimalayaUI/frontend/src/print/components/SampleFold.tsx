import { useState, useRef, useEffect, type JSX } from "react";
import type { LoadSample } from "../../api";
import { Checkbox } from "../ui/Checkbox";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";
import { ExposureLeaf } from "./ExposureLeaf";

export interface SampleFoldProps {
  sample: LoadSample;
  open: boolean;
  selected: boolean;
  onToggleOpen: (sampleId: number) => void;
  onToggleSelect: (sampleId: number) => void;
  /** Called with the sample id and trimmed new name when the rename is committed. */
  onRename: (sampleId: number, newName: string) => void;
  onSplit: (sampleId: number) => void;
  /** Merge this sample (loser) into the flagged partner (survivor). */
  onMerge: (loserId: number, survivorId: number) => void;
  onDismissFlag: (sampleId: number) => void;
  onMoveExposure: (sampleId: number, exposureId: number, anchorEl: HTMLElement) => void;
  /** Per-exposure thumbnail URL resolver (page supplies; Thumbnail takes a src). */
  thumbSrcFor: (exposureId: number) => string | null;
  className?: string;
}

export function SampleFold(props: SampleFoldProps): JSX.Element {
  const { sample: s, open, selected } = props;
  const flag = s.flag;
  const flagChip =
    flag?.kind === "split" ? "check split" :
    flag?.kind === "merge" ? "possible reshoot" : null;

  // Inline rename state: null = not editing; string = draft value being edited.
  const [renameDraft, setRenameDraft] = useState<string | null>(null);
  const renameInputRef = useRef<HTMLInputElement>(null);
  const isRenaming = renameDraft !== null;

  // Focus the input whenever rename mode opens.
  useEffect(() => {
    if (isRenaming) renameInputRef.current?.select();
  }, [isRenaming]);

  const activateRename = () => setRenameDraft(s.name);

  const commitRename = () => {
    if (renameDraft === null) return;
    const trimmed = renameDraft.trim();
    setRenameDraft(null);
    if (trimmed && trimmed !== s.name) {
      props.onRename(s.sample_id, trimmed);
    }
  };

  const cancelRename = () => setRenameDraft(null);

  return (
    <div
      data-testid="sample-fold"
      className={[
        "rounded-sm border border-hair bg-plate",
        flag ? "border-warning" : "",
        selected ? "ring-1 ring-accent" : "",
        props.className ?? "",
      ].filter(Boolean).join(" ")}
    >
      <div className="flex items-center gap-3 px-3.5 py-2.5">
        <Checkbox
          checked={selected}
          onChange={() => props.onToggleSelect(s.sample_id)}
          aria-label={`Select ${s.name}`}
        />
        {isRenaming ? (
          <div className="min-w-0 flex-1">
            <Input
              variant="title"
              testId="sample-rename-input"
              value={renameDraft}
              onValueChange={setRenameDraft}
              inputRef={renameInputRef}
              aria-label="Sample name"
              onBlur={commitRename}
              onKeyDown={(e) => {
                if (e.key === "Enter") {
                  commitRename();
                  (e.target as HTMLElement).blur();
                } else if (e.key === "Escape") {
                  cancelRename();
                  (e.target as HTMLElement).blur();
                }
              }}
            />
          </div>
        ) : (
          <button
            type="button"
            className="min-w-0 flex-1 text-left"
            onClick={() => props.onToggleOpen(s.sample_id)}
          >
            <span className="text-headline text-ink">{s.name}</span>
            <span className="ml-2 font-mono text-xs text-ink-faint">
              {s.exposures.length} exposures · {s.exposures[0]?.timestamp ?? "--"}
            </span>
            {flagChip ? (
              <span className="ml-2 text-xs font-bold uppercase text-warning">{flagChip}</span>
            ) : null}
          </button>
        )}
        <div className="flex items-center gap-1.5">
          <Button variant="outline" onClick={activateRename}>Rename</Button>
          <Button variant="outline" onClick={() => props.onSplit(s.sample_id)}>Split...</Button>
        </div>
      </div>

      {open ? (
        <div className="border-t border-hair p-1">
          {flag?.kind === "merge" ? (
            <div
              data-testid="merge-prompt"
              className="m-2.5 rounded-sm border border-warning bg-paper-sunk px-3.5 py-2.5"
            >
              <div className="text-xs text-ink-soft">
                Filename matches <b className="text-ink">{flag.merge_with_label}</b>. This looks like a reshoot in a later load.
              </div>
              <div className="mt-2 flex gap-2">
                <Button
                  variant="outline"
                  onClick={() => props.onMerge(s.sample_id, flag.merge_with_sample_id)}
                >
                  Merge into that sample
                </Button>
                <Button
                  variant="ghost"
                  onClick={() => props.onDismissFlag(s.sample_id)}
                >
                  Keep separate
                </Button>
              </div>
            </div>
          ) : null}

          {s.exposures.map((e, i) => (
            <div key={e.id}>
              {flag?.kind === "split" && i === flag.split_at_index ? (
                <div
                  data-testid="split-divider"
                  className="mx-2.5 my-1 flex items-center gap-3 rounded-sm border border-warning bg-paper-sunk px-3 py-1.5"
                >
                  <span className="flex-1 font-mono text-xs text-ink-soft">
                    position jumps{" "}
                    <b className="text-warning">
                      {flag.jump_from.toFixed(1)} {"→"} {flag.jump_to.toFixed(1)} mm
                    </b>
                    {" "}here -- likely two samples
                  </span>
                  <Button variant="outline" onClick={() => props.onSplit(s.sample_id)}>
                    Split here
                  </Button>
                </div>
              ) : null}
              <ExposureLeaf
                exposure={e}
                thumbSrc={props.thumbSrcFor(e.id)}
                onMove={(eid, anchorEl) => props.onMoveExposure(s.sample_id, eid, anchorEl)}
              />
            </div>
          ))}
        </div>
      ) : null}
    </div>
  );
}
