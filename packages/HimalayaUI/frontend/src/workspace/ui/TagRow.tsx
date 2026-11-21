import type { SampleLabel } from "../../api";
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
        return (
          <span
            key={t.id}
            data-tag={t.label}
            data-tag-state={t.state}
            className={
              "inline-flex items-center gap-1 rounded-sm border bg-plate py-0.5 pl-2 pr-1 text-caption " +
              (suggested ? "border-dashed border-ink-faint text-ink-soft" : "border-hair text-ink")
            }
          >
            {suggested && editable ? (
              <button type="button" onClick={() => onAdd(t.label)} title="Suggested. Click to confirm" aria-label={`Confirm tag ${t.label}`}>
                {t.label}?
              </button>
            ) : (
              <span>{suggested ? `${t.label}?` : t.label}</span>
            )}
            {editable && (
              <button type="button" aria-label={`Remove tag ${t.label}`} onClick={() => onRemove(t.label)} className="px-1 text-ink-faint hover:text-ink">
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
          className="inline-flex items-center rounded-sm border border-dashed border-hair-strong px-2 py-0.5 text-caption text-ink-soft hover:text-ink"
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
