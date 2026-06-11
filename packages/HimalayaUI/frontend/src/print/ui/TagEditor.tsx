import { useId, useState } from "react";
import type { KeyboardEvent } from "react";
import { Input } from "./Input";
import { Button } from "./Button";
import { Chip } from "./Chip";
import type { Tag } from "./tag";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

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
    <div
      data-testid="tag-editor"
      className={cx("inline-flex items-center gap-1.5", className)}
    >
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
  );
}
