import { useEffect, useState } from "react";
import type { Sample } from "../api";

interface ExposureSummary {
  total: number;
  accepted: number;
  rejected: number;
}

interface Props {
  sample: Sample;
  experimentName?: string | undefined;
  exposureSummary: ExposureSummary;
  onUpdateSample: (patch: { display_name?: string; notes?: string }) => void;
  onAddTag: (key: string, value: string) => void;
  onRemoveTag: (tagId: number) => void;
}

export function SampleMetadataCard({
  sample,
  experimentName,
  exposureSummary,
  onUpdateSample,
  onAddTag,
  onRemoveTag,
}: Props): JSX.Element {
  const [name,  setName]  = useState(sample.display_name  ?? "");
  const [notes, setNotes] = useState(sample.notes ?? "");
  const [nameFocused,  setNameFocused]  = useState(false);
  const [notesFocused, setNotesFocused] = useState(false);
  const [newTagKey, setNewTagKey] = useState("");
  const [newTagVal, setNewTagVal] = useState("");
  const [addingTag, setAddingTag] = useState(false);

  // Re-sync editable fields when the active sample (or its server-side
  // values) change. useState initializers only fire on first mount, so
  // without this, switching samples leaves the previous sample's values
  // in the inputs.
  //
  // Focus-gate (mirrors QNumInput in PlotCard.tsx): skip the sync while
  // the user is mid-edit, so a background refetch or another user's save
  // doesn't clobber an in-progress draft. Sample-switch still works
  // because focus moves out before the new sample renders.
  useEffect(() => {
    if (!nameFocused)  setName(sample.display_name ?? "");
    if (!notesFocused) setNotes(sample.notes ?? "");
  }, [sample.id, sample.display_name, sample.notes, nameFocused, notesFocused]);

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
    <div className="flex flex-col gap-3 p-3 overflow-y-auto">
      {/* Breadcrumb + name — leads the card */}
      <div className="flex flex-col gap-1">
        {/* experimentName · stable identifier breadcrumb */}
        {(experimentName || sample.name) && (
          <p className="text-caption truncate">
            {experimentName && (
              <span>{experimentName}</span>
            )}
            {experimentName && sample.name && (
              <span className="mx-1 opacity-40">·</span>
            )}
            {sample.name && (
              <span className="font-medium text-fg-dim">{sample.name}</span>
            )}
          </p>
        )}
        {/* data-testid is historical (issue #88); the field now edits display_name. */}
        <input
          data-testid="sample-name-input"
          className="w-full bg-transparent border-0 outline-none px-0 text-title placeholder:text-fg-muted"
          placeholder="Sample name"
          value={name}
          onChange={(e) => setName(e.target.value)}
          onFocus={() => setNameFocused(true)}
          onBlur={() => {
            setNameFocused(false);
            onUpdateSample({ display_name: name });
          }}
        />
      </div>

      <p className="text-caption">
        {exposureSummary.total} exposures
        {" · "}
        {exposureSummary.accepted} accepted
        {" · "}
        {exposureSummary.rejected} rejected
      </p>

      {/* Notes */}
      <div className="flex flex-col gap-0.5">
        <label className="text-label">
          Notes
        </label>
        <textarea
          data-testid="sample-notes-textarea"
          rows={2}
          className="w-full bg-bg border border-border rounded px-2 py-1 text-sm text-fg resize-none"
          value={notes}
          onChange={(e) => setNotes(e.target.value)}
          onFocus={() => setNotesFocused(true)}
          onBlur={() => {
            setNotesFocused(false);
            onUpdateSample({ notes });
          }}
        />
      </div>

      {/* Tags */}
      <div className="flex flex-col gap-1.5">
        <label className="text-label">
          Tags
        </label>
        <div className="flex flex-wrap gap-1">
          {sample.tags.map((tag) => (
            <span
              key={tag.id}
              className="inline-flex items-center gap-1 text-xs px-2 py-0.5 rounded-full
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
              className="text-xs text-fg-muted hover:text-fg px-2 py-0.5 rounded-full
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
