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

export function TagRow({ tags, vocabulary, editable, adding, onOpenAdd, onCloseAdd, onAdd, onRemove }: Props): JSX.Element | null {
  if (tags.length === 0 && !editable) return null;
  const last = tags.filter((t) => t.state === "confirmed").at(-1);
  return (
    <div className="relative flex flex-wrap items-center gap-1.5 px-3.5 pb-1" data-testid="workspace-tags">
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
            className={`inline-flex items-center gap-0.5 rounded-full border py-0.5 pl-2.5 pr-1 text-caption ${tone} ${fill}`}
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
      {editable && (
        <button
          type="button"
          onClick={onOpenAdd}
          className="inline-flex items-center rounded-full border border-dashed border-hair-strong px-2.5 py-0.5 text-caption text-ink-soft hover:text-ink"
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
