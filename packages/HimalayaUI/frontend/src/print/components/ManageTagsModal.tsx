import { useEffect, useId, useRef, useState } from "react";
import type { RefObject } from "react";
import { ModalShell, IconButton, Button } from "../ui";
import type { TagSuggestOption } from "../ui";
import { TagSuggest } from "../ui";
import { ModalHead } from "./ModalHead";
import { ModalFooter } from "./ModalFooter";
import type { LoupeTag } from "../pages/loupeAdapters";
import { announce } from "../../lib/announce";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ManageTagsModalProps {
  open: boolean;
  sampleName: string;
  tags: LoupeTag[];
  keyOptions: TagSuggestOption[];
  valueOptionsFor: (key: string) => TagSuggestOption[];
  onEdit: (tagId: number, key: string, value: string) => void;
  onAdd: (key: string, value: string) => void;
  onRemove: (tagId: number) => void;
  onClose: () => void;
  /** Ref to the element that opened this modal — focus is restored to it on close. */
  triggerRef?: RefObject<HTMLButtonElement>;
}

interface RowDraft {
  key: string;
  value: string;
}

/**
 * ManageTagsModal — per-sample tag editor opened from the "Manage" affordance
 * in LoupeSidePanel. Design-guard SCANNED (lives in src/print/components/) →
 * token-only classes, mirrors CustomIndexModal structure.
 *
 * Each existing tag renders as an editable row: key TagSuggest + value
 * TagSuggest + dismiss IconButton. A trailing add-row lets the user type a
 * new key+value pair. Inline `role="alert"` blocks duplicate-key commits
 * (single-valued-key rule). Focus is restored to the `triggerRef` on close
 * (ModalShell does not do this). Announces add/edit/remove via lib/announce.
 */
export function ManageTagsModal({
  open,
  sampleName,
  tags,
  keyOptions,
  valueOptionsFor,
  onEdit,
  onAdd,
  onRemove,
  onClose,
  triggerRef,
}: ManageTagsModalProps): JSX.Element | null {
  const titleId = useId();
  const addErrorId = useId();
  const editErrorIdPrefix = useId();

  // Per-row draft state (key + value text mirrors) — keyed by tag id
  const [drafts, setDrafts] = useState<Record<number, RowDraft>>({});
  // The tag id whose key-edit collided with another row's key (inline reject).
  const [editKeyDupId, setEditKeyDupId] = useState<number | null>(null);
  // Add-row state
  const [addKey, setAddKey] = useState("");
  const [addValue, setAddValue] = useState("");
  const [addKeyDup, setAddKeyDup] = useState(false);

  // Seed drafts when tags change (open or new tags arrive)
  useEffect(() => {
    if (!open) return;
    setDrafts((prev) => {
      const next: Record<number, RowDraft> = {};
      for (const t of tags) {
        next[t.id] = prev[t.id] ?? { key: t.key, value: t.value ?? "" };
      }
      return next;
    });
  }, [open, tags]);

  // Reset add-row on open
  useEffect(() => {
    if (!open) {
      setAddKey("");
      setAddValue("");
      setAddKeyDup(false);
      setEditKeyDupId(null);
      setDrafts({});
    }
  }, [open]);

  // Restore focus to the trigger ONLY on a genuine open→close transition.
  // Without the wasOpenRef guard this effect also fires on initial mount (the
  // modal mounts with open=false on every loupe load) and would steal focus to
  // the "Manage" button on page load — a real focus-steal bug.
  const wasOpenRef = useRef(open);
  useEffect(() => {
    if (wasOpenRef.current && !open) {
      triggerRef?.current?.focus();
    }
    wasOpenRef.current = open;
  }, [open, triggerRef]);

  if (!open) return null;

  const existingKeys = tags.map((t) => t.key);

  const setDraft = (id: number, field: "key" | "value", v: string): void => {
    setDrafts((prev) => ({
      ...prev,
      [id]: { ...(prev[id] ?? { key: "", value: "" }), [field]: v },
    }));
    // Editing the key of the row that's flagged clears its collision error.
    if (field === "key" && editKeyDupId === id) setEditKeyDupId(null);
  };

  const commitEdit = (tag: LoupeTag, field: "key" | "value", v: string): void => {
    const draft = drafts[tag.id] ?? { key: tag.key, value: tag.value ?? "" };
    const next = field === "key" ? { ...draft, key: v } : { ...draft, value: v };
    const nextKey = next.key.trim();
    // Single-valued-key rule on the EDIT path: reject inline if the new key
    // collides with ANOTHER row's key on the same sample (excluding this row).
    // The backend 409s on this too, but the spec requires inline rejection.
    if (field === "key") {
      const collides = tags.some((t) => t.id !== tag.id && t.key === nextKey);
      if (collides) {
        setDrafts((prev) => ({ ...prev, [tag.id]: next }));
        setEditKeyDupId(tag.id);
        return;
      }
    }
    // Update local draft immediately (TagSuggest fires onCommit)
    setDrafts((prev) => ({ ...prev, [tag.id]: next }));
    if (editKeyDupId === tag.id) setEditKeyDupId(null);
    onEdit(tag.id, nextKey, next.value.trim());
    announce(`Tag updated: ${next.key}${next.value ? ` ${next.value}` : ""}`);
  };

  const handleRemove = (tag: LoupeTag): void => {
    onRemove(tag.id);
    announce(`Tag removed: ${tag.key}${tag.value ? ` ${tag.value}` : ""}`);
  };

  const handleAdd = (): void => {
    const k = addKey.trim();
    const v = addValue.trim();
    if (!k) return;
    if (existingKeys.includes(k)) {
      setAddKeyDup(true);
      return;
    }
    onAdd(k, v);
    announce(`Tag added: ${k}${v ? ` ${v}` : ""}`);
    setAddKey("");
    setAddValue("");
    setAddKeyDup(false);
  };

  const clearAddDupError = (newKey: string): void => {
    setAddKey(newKey);
    if (addKeyDup) setAddKeyDup(false);
  };

  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="md"
      testId="manage-tags-modal"
      aria-labelledby={titleId}
    >
      <ModalHead
        kicker="Sample tags"
        title={sampleName}
        titleId={titleId}
        onClose={onClose}
      />

      {/* Tag rows — shared 3-column grid: [key][value][action].
          All rows (edit + add) sit in the same grid so the action column
          auto-sizes to the widest child (the "Add" button) and × aligns
          with it above. `relative z-10` creates a stacking context above the
          footer so TagSuggest dropdowns (z-30 within this context) render
          on top of the footer copy. */}
      <div className="px-5 pt-4 pb-2 relative z-10">
        <div className="grid grid-cols-[1fr_1fr_auto] gap-x-2 gap-y-3">
          {tags.map((tag) => {
            const draft = drafts[tag.id] ?? { key: tag.key, value: tag.value ?? "" };
            const keyDup = editKeyDupId === tag.id;
            const rowErrorId = `${editErrorIdPrefix}-${tag.id}`;
            return (
              <div key={tag.id} className="contents">
                <div
                  data-testid="manage-tag-row"
                  className="col-span-3 grid grid-cols-subgrid items-center"
                >
                  <TagSuggest
                    label={`Key for tag ${tag.id}`}
                    mode="key"
                    value={draft.key}
                    options={keyOptions}
                    onValueChange={(v) => setDraft(tag.id, "key", v)}
                    onCommit={(v) => commitEdit(tag, "key", v)}
                    invalid={keyDup}
                    {...(keyDup ? { "aria-describedby": rowErrorId } : {})}
                    className="min-w-0 w-full"
                  />
                  <TagSuggest
                    label={`Value for tag ${tag.id}`}
                    mode="value"
                    value={draft.value}
                    options={valueOptionsFor(draft.key)}
                    onValueChange={(v) => setDraft(tag.id, "value", v)}
                    onCommit={(v) => commitEdit(tag, "value", v)}
                    className="min-w-0 w-full"
                  />
                  <IconButton
                    dismiss
                    label={`Remove ${tag.key}${tag.value ? ` ${tag.value}` : ""}`}
                    tone="ghost"
                    onClick={() => handleRemove(tag)}
                  />
                </div>
                {keyDup && (
                  <span
                    id={rowErrorId}
                    role="alert"
                    data-testid="manage-tag-edit-dup-error"
                    className="col-span-3 text-caption text-error"
                  >
                    This sample already has a tag with that key.
                  </span>
                )}
              </div>
            );
          })}

          {/* Add-a-tag row — aligns with the edit rows above */}
          <div
            data-testid="manage-tag-add-row"
            className="col-span-3 grid grid-cols-subgrid items-center"
          >
            <TagSuggest
              label="New tag key"
              mode="key"
              value={addKey}
              options={keyOptions}
              onValueChange={clearAddDupError}
              onCommit={(v) => {
                clearAddDupError(v);
              }}
              invalid={addKeyDup}
              {...(addKeyDup ? { "aria-describedby": addErrorId } : {})}
              className="min-w-0 w-full"
            />
            <TagSuggest
              label="New tag value"
              mode="value"
              value={addValue}
              options={valueOptionsFor(addKey)}
              onValueChange={setAddValue}
              onCommit={(v) => setAddValue(v)}
              className="min-w-0 w-full"
            />
            <Button variant="solid" onClick={handleAdd}>
              Add
            </Button>
          </div>
          {addKeyDup && (
            <span
              id={addErrorId}
              role="alert"
              data-testid="manage-tag-dup-error"
              className={cx("col-span-3 text-caption text-error")}
            >
              This sample already has a tag with that key.
            </span>
          )}
        </div>
      </div>

      <ModalFooter
        note="Tags apply to the whole sample. Saved as you go."
        actions={
          <Button variant="outline" onClick={onClose}>
            Done
          </Button>
        }
      />
    </ModalShell>
  );
}
