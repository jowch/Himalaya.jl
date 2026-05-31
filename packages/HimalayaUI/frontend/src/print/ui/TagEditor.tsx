import { useState } from "react";
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
 *  C — focus/hover/disabled flow through the composed Input/Button primitives.
 *  F — appearance lives in those primitives; this component contributes
 *  placement-only layout (`inline-flex`, gaps). */
export function TagEditor({
  onCommit,
  onCancel,
  knownKeys,
  className = "",
}: TagEditorProps): JSX.Element {
  const [key, setKey] = useState("");
  const [value, setValue] = useState("");

  const commit = (): void => {
    const k = key.trim();
    if (!k) return;
    const v = value.trim();
    onCommit({ key: k, ...(v ? { value: v } : {}) });
    setKey("");
    setValue("");
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
        onValueChange={setKey}
        placeholder="key"
        inputSize="sm"
        onKeyDown={onKeyDown}
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
      {knownKeys && knownKeys.length > 0 ? (
        <span className="inline-flex items-center gap-1">
          {knownKeys.map((kk) => (
            <Chip key={kk} variant="add" onClick={() => setKey(kk)}>
              {kk}
            </Chip>
          ))}
        </span>
      ) : null}
    </div>
  );
}
