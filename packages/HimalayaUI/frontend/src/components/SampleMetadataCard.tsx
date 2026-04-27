import { useState } from "react";
import type { Sample } from "../api";

interface ExposureSummary {
  total: number;
  accepted: number;
  rejected: number;
}

interface Props {
  sample: Sample;
  exposureSummary: ExposureSummary;
  onUpdateSample: (patch: { name?: string; notes?: string }) => void;
  onAddTag: (key: string, value: string) => void;
  onRemoveTag: (tagId: number) => void;
}

export function SampleMetadataCard({
  sample,
  exposureSummary,
  onUpdateSample,
  onAddTag,
  onRemoveTag,
}: Props): JSX.Element {
  const [name,  setName]  = useState(sample.name  ?? "");
  const [notes, setNotes] = useState(sample.notes ?? "");
  const [newTagKey, setNewTagKey] = useState("");
  const [newTagVal, setNewTagVal] = useState("");
  const [addingTag, setAddingTag] = useState(false);

  function handleAddTag() {
    const k = newTagKey.trim();
    const v = newTagVal.trim();
    if (k && v) {
      onAddTag(k, v);
      setNewTagKey("");
      setNewTagVal("");
      setAddingTag(false);
    }
  }

  return (
    <div className="card flex flex-col gap-3 p-3 overflow-y-auto">
      {/* Name — leads the card */}
      <div className="flex flex-col gap-0.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">
          Name
        </label>
        <input
          className="w-full bg-bg border border-border rounded px-2 py-1 text-[15px] font-semibold text-fg"
          value={name}
          onChange={(e) => setName(e.target.value)}
          onBlur={() => onUpdateSample({ name })}
        />
      </div>

      <p className="text-[10px] text-fg-muted">
        {exposureSummary.total} exposures
        {" · "}
        {exposureSummary.accepted} accepted
        {" · "}
        {exposureSummary.rejected} rejected
      </p>

      {/* Notes */}
      <div className="flex flex-col gap-0.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">
          Notes
        </label>
        <textarea
          rows={2}
          className="w-full bg-bg border border-border rounded px-2 py-1 text-sm text-fg resize-none"
          value={notes}
          onChange={(e) => setNotes(e.target.value)}
          onBlur={() => onUpdateSample({ notes })}
        />
      </div>

      {/* Label (read-only) */}
      {sample.label && (
        <div className="flex flex-col gap-0.5">
          <label className="text-[10px] text-fg-muted uppercase tracking-wide">
            Label
          </label>
          <span className="text-xs text-fg-muted bg-bg-subtle px-2 py-1 rounded">
            {sample.label}
          </span>
        </div>
      )}

      {/* Tags */}
      <div className="flex flex-col gap-1.5">
        <label className="text-[10px] text-fg-muted uppercase tracking-wide">
          Tags
        </label>
        <div className="flex flex-wrap gap-1">
          {sample.tags.map((tag) => (
            <span
              key={tag.id}
              className="inline-flex items-center gap-1 text-[10px] px-2 py-0.5 rounded-full
                         bg-bg-subtle border border-border text-fg-muted"
            >
              {tag.key}: {tag.value}
              {tag.source !== "manifest" && (
                <button
                  onClick={() => onRemoveTag(tag.id)}
                  className="text-fg-muted hover:text-error leading-none"
                  aria-label={`Remove ${tag.key} tag`}
                >
                  ×
                </button>
              )}
            </span>
          ))}
          {!addingTag && (
            <button
              onClick={() => setAddingTag(true)}
              className="text-[10px] text-fg-muted hover:text-fg px-2 py-0.5 rounded-full
                         border border-dashed border-border"
            >
              + tag
            </button>
          )}
        </div>
        {addingTag && (
          <div className="flex gap-1">
            <input
              className="flex-1 text-xs bg-bg border border-border rounded px-1.5 py-0.5"
              placeholder="key"
              value={newTagKey}
              onChange={(e) => setNewTagKey(e.target.value)}
              autoFocus
            />
            <input
              className="flex-1 text-xs bg-bg border border-border rounded px-1.5 py-0.5"
              placeholder="value"
              value={newTagVal}
              onChange={(e) => setNewTagVal(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && handleAddTag()}
            />
            <button
              onClick={handleAddTag}
              className="text-xs px-2 border border-accent text-accent rounded"
            >
              Add
            </button>
            <button
              onClick={() => setAddingTag(false)}
              className="text-xs px-2 text-fg-muted"
            >
              ×
            </button>
          </div>
        )}
      </div>
    </div>
  );
}
