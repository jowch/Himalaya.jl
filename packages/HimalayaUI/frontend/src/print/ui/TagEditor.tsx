import { useId, useState } from "react";
import type { KeyboardEvent } from "react";
import { Input } from "./Input";
import { Button } from "./Button";
import { Chip } from "./Chip";
import type { Tag } from "./tag";
import { cx } from "../../lib/cx";


export interface TagEditorProps {
  onCommit: (tag: Tag) => void;
  onCancel?: () => void;
  /** Optional assist: known tag keys offered as one-click suggestions to fill
   *  the key field. Free-text key entry is always available regardless. */
  knownKeys?: string[];
  /** Keys already on this sample — adding one of these is rejected inline
   *  (single-valued-key rule). When undefined the check is skipped (backward-
   *  compatible: all existing consumers that don't pass this still work). */
  existingKeys?: string[];
  /** PLACEMENT-ONLY: margin / alignment. No appearance utils. */
  className?: string;
}

/** Inline key(+optional value) tag entry (DESIGN.md tag model — backend tags are
 *  key-value, some key-only). Composes the base `Input` twice plus a solid
 *  commit `Button`. Enter (in either field) or clicking "Add" commits; Escape
 *  cancels. An empty value yields a key-only tag via conditional spread (no
 *  `value` key at all). `knownKeys`, when present, render as `Chip variant="add"`
 *  suggestions that fill the key field.
 *
 *  Committing with an EMPTY key is rejected honestly, not silently: the key
 *  field flips `aria-invalid` (via Input's `invalid` border) and a quiet
 *  inline error appears, tied to the field by `aria-describedby`. The error
 *  clears as soon as a key is typed. Validation lives here because the
 *  primitive owns commit.
 *
 *  C — focus/hover/disabled flow through the composed Input/Button primitives.
 *  F — appearance lives in those primitives; this component contributes
 *  placement-only layout (`inline-flex`, gaps). */
export function TagEditor({
  onCommit,
  onCancel,
  knownKeys,
  existingKeys,
  className = "",
}: TagEditorProps): JSX.Element {
  const [key, setKey] = useState("");
  const [value, setValue] = useState("");
  const [keyMissing, setKeyMissing] = useState(false);
  const [keyDuplicate, setKeyDuplicate] = useState(false);
  const errorId = useId();

  const commit = (): void => {
    const k = key.trim();
    if (!k) {
      // Honest rejection: an empty-key Add is invalid input, never a silent no-op.
      setKeyMissing(true);
      setKeyDuplicate(false);
      return;
    }
    if (existingKeys?.includes(k)) {
      // Single-valued-key rule: reject adding a key the sample already has.
      setKeyDuplicate(true);
      setKeyMissing(false);
      return;
    }
    const v = value.trim();
    onCommit({ key: k, ...(v ? { value: v } : {}) });
    setKey("");
    setValue("");
    setKeyMissing(false);
    setKeyDuplicate(false);
  };

  const onKeyChange = (v: string): void => {
    setKey(v);
    if ((keyMissing || keyDuplicate) && v.trim() !== "") {
      setKeyMissing(false);
      setKeyDuplicate(false);
    }
  };

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.key === "Enter") {
      e.preventDefault();
      commit();
    } else if (e.key === "Escape") {
      e.preventDefault();
      onCancel?.();
    }
  };

  return (
    // TAG-ERR-LAYOUT: a column — the controls wrap on one row, the validation
    // error renders as its own full-width line BELOW them (in the error token).
    // Previously the error sat inside the inline controls row and got squeezed
    // into the leftover ~40px, stacking one word per line at the 286px rail.
    <div
      data-testid="tag-editor"
      // data-keys-trap: while this inline editor is open, page-level single-key
      // shortcuts (loupe X/R, sheet X) must not fire through its buttons/chips
      // — suppressGlobalKeys treats the subtree as a keyboard island (LO-POPKEY).
      data-keys-trap=""
      className={cx("flex flex-col gap-1 items-start", className)}
    >
      <div className="flex flex-wrap items-center gap-1.5">
        <Input
          value={key}
          onValueChange={onKeyChange}
          placeholder="key"
          inputSize="sm"
          onKeyDown={onKeyDown}
          invalid={keyMissing || keyDuplicate}
          {...((keyMissing || keyDuplicate) ? { "aria-describedby": errorId } : {})}
        />
        <Input
          value={value}
          onValueChange={setValue}
          placeholder="value (optional)"
          inputSize="sm"
          onKeyDown={onKeyDown}
        />
        <Button variant="solid" onClick={commit}>
          Add
        </Button>
        {onCancel ? (
          // Visible escape hatch out of the inline add (Escape also cancels, but
          // a pointer user needs a target). Quiet `ghost` so it sits beside the
          // solid Add without competing for the eye.
          <Button variant="ghost" onClick={onCancel}>
            Cancel
          </Button>
        ) : null}
        {knownKeys && knownKeys.length > 0 ? (
          <span className="inline-flex items-center gap-1">
            {knownKeys.map((kk) => (
              <Chip key={kk} variant="add" onClick={() => onKeyChange(kk)}>
                {kk}
              </Chip>
            ))}
          </span>
        ) : null}
      </div>
      {keyMissing ? (
        <span
          id={errorId}
          role="alert"
          data-testid="tag-editor-error"
          className="text-caption text-error"
        >
          Enter a key first.
        </span>
      ) : keyDuplicate ? (
        <span
          id={errorId}
          role="alert"
          data-testid="tag-editor-error"
          className="text-caption text-error"
        >
          This sample already has a tag with that key.
        </span>
      ) : null}
    </div>
  );
}
