import type { SampleLabel } from "../../api";
import { displayTag, standardTag } from "../../lib/labels";
import { TagPicker } from "./TagPicker";

interface Props {
  tags: readonly SampleLabel[];
  vocabulary: readonly string[];
  editable: boolean;
  adding: boolean;
  onOpenAdd: () => void;
  onCloseAdd: () => void;
  onAdd: (label: string) => void;
  onRemove: (label: string) => void;
}

export function TagRow({ tags, vocabulary, editable, adding, onOpenAdd, onCloseAdd, onAdd, onRemove }: Props): JSX.Element {
  const last = tags.filter((t) => t.state === "confirmed").at(-1);
  return (
    <div className="relative flex h-7 flex-nowrap items-center gap-1.5 px-3.5" data-testid="workspace-tags">
      <div className="flex min-w-0 flex-nowrap items-center gap-1.5 overflow-hidden" data-testid="workspace-tag-chips">
        {tags.map((t) => {
          const suggested = t.state === "suggested";
          const standard = standardTag(t.label) !== undefined;
          const name = displayTag(t.label);
          const tone = standard ? "text-tag-std-ink" : "text-tag-ink";
          const fill = suggested
            ? "border-dashed border-current bg-transparent"
            : standard
              ? "border-transparent bg-tag-std-wash"
              : "border-transparent bg-tag-wash";
          return (
            <span
              key={t.id}
              data-tag={t.label}
              data-tag-state={t.state}
              data-tag-standard={standard ? "true" : undefined}
              className={`inline-flex shrink-0 items-center gap-0.5 whitespace-nowrap rounded-full border pb-0 pl-2.5 pr-1 pt-1 text-caption ${tone} ${fill}`}
            >
              {suggested && editable ? (
                <button type="button" onClick={() => onAdd(t.label)} title="Suggested. Click to confirm" aria-label={`Confirm tag ${name}`}>
                  {name}
                </button>
              ) : (
                <span>{name}</span>
              )}
              {editable && (
                <button
                  type="button"
                  aria-label={`Remove tag ${name}`}
                  onClick={() => onRemove(t.label)}
                  className="rounded-full px-1 opacity-60 hover:opacity-100"
                >
                  ×
                </button>
              )}
            </span>
          );
        })}
      </div>
      {editable && (
        <button
          type="button"
          onClick={onOpenAdd}
          className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-dashed border-hair-strong px-2.5 pb-0 pt-1 text-caption text-ink-soft hover:text-ink"
        >
          + Tag
        </button>
      )}
      {adding && editable && (
        <TagPicker
          vocabulary={vocabulary}
          onSample={tags.filter((t) => t.state === "confirmed").map((t) => t.label)}
          onPick={onAdd}
          onRemoveLast={() => {
            if (last) onRemove(last.label);
          }}
          onClose={onCloseAdd}
        />
      )}
    </div>
  );
}
