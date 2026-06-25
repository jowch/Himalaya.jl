import { type JSX } from "react";
import { useInlineEdit } from "../../hooks/useInlineEdit";
import type { LoadSample } from "../../api";
import { Checkbox } from "../ui/Checkbox";
import { Button } from "../ui/Button";
import { IconButton } from "../ui/IconButton";
import { Input } from "../ui/Input";
import { ExposureLeaf } from "./ExposureLeaf";

export interface SampleFoldProps {
  sample: LoadSample;
  open: boolean;
  selected: boolean;
  /** Keyboard cursor lands here (↑/↓ nav). Distinct from `selected`: an outline
   *  (focus look) vs the selected accent ring + checked box. */
  cursored?: boolean;
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
  // 0-based split boundary (start of the split-off tail), clamped to [1, len-1]
  // — the SAME boundary GroupingReviewPage.handleSplit uses, so the divider
  // always sits exactly where the split would land (even for an out-of-range
  // flag the backend shouldn't emit). -1 when not a split flag.
  const splitBoundary =
    flag?.kind === "split"
      ? Math.min(Math.max(flag.split_at_index - 1, 1), s.exposures.length - 1)
      : -1;

  // Inline rename: the hook focuses+selects on open, trims, and skips a no-op.
  // A blank rename is ignored here (never blank a sample's name).
  const {
    editingKey,
    draft: renameDraft,
    setDraft: setRenameDraft,
    inputRef: renameInputRef,
    begin,
    commit: commitRename,
    cancel: cancelRename,
  } = useInlineEdit<true>((_key, trimmed) => {
    if (trimmed) props.onRename(s.sample_id, trimmed);
  });
  const isRenaming = editingKey !== null;

  return (
    <div
      data-testid="sample-fold"
      data-cursored={props.cursored ? "true" : undefined}
      className={[
        // scroll-mb-14 (56px): reserves clearance below for the page's
        // scrollIntoView so a downward-cursored row aligns ABOVE the fixed
        // grouping-footer Dock (~47px), not behind it. Mirrors SampleTableRow.
        "rounded-sm border border-hair bg-plate scroll-mb-14",
        flag ? "border-warning" : "",
        selected ? "ring-1 ring-accent" : "",
        // Cursor = an offset outline (focus look), independent of the selected
        // ring (box-shadow) so the two can show at once without clobbering.
        props.cursored ? "outline outline-2 outline-offset-2 outline-accent" : "",
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
              inputSize="sm"
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
          <div className="min-w-0 flex-1 flex items-center gap-1.5">
            <button
              type="button"
              className="min-w-0 text-left"
              onClick={() => props.onToggleOpen(s.sample_id)}
            >
              <span className="text-headline text-ink">{s.name}</span>
              <span className="ml-2 font-mono text-xs text-ink-faint">
                {s.exposures.length} exposures · {s.exposures[0]?.timestamp ?? "--"}
              </span>
            </button>
            {/* Pencil sits next to the name (the funnel rename idiom), not as a
                far-right button. */}
            <IconButton label="Rename sample" tone="ghost" onClick={() => begin(true, s.name)}>
              {"✎"}
            </IconButton>
          </div>
        )}
        {/* Flag chip + an explicit dismiss (×), discoverable WITHOUT expanding,
            for BOTH split and merge kinds. The backend dismiss route handles
            both; split flags previously had no dismiss affordance at all. */}
        {flag && flagChip ? (
          <div className="flex items-center gap-1">
            <span className="text-xs font-bold uppercase text-warning">{flagChip}</span>
            <IconButton
              label="Dismiss flag"
              tone="ghost"
              dismiss
              onClick={() => props.onDismissFlag(s.sample_id)}
            />
          </div>
        ) : null}
        <div className="flex items-center gap-1.5">
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
              {/* The divider sits at the clamped 0-based split boundary — the
                  start of the split-off tail (== handleSplit's slice point). */}
              {flag?.kind === "split" && i === splitBoundary ? (
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
