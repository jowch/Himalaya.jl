import { useEffect, useRef, useState } from "react";

interface Props {
  value: string;
  onCommit: (newValue: string) => void;
  /** Optional placeholder when value is empty */
  placeholder?: string;
  /** className for the rest-state element */
  className?: string;
  /** className for the input */
  inputClassName?: string;
  /** Test selector */
  testId?: string;
  /**
   * When true, renders as plain text — no click handler, no hover affordance,
   * no `data-interactable`. Used by review-mode `CompareTitleStrip` to render
   * the title as read-only without silently dropping clicks. Default `false`.
   */
  readOnly?: boolean;
}

/**
 * Click text → input. Enter commits, Esc cancels (restores prior value),
 * blur commits. Carries `data-interactable="edit"` for the
 * visual-language vocabulary (spec §4).
 *
 * When `readOnly` is true, behaves as a static `<span>` with no interactive
 * affordance. Caller can flip this based on `useCompareMode().kind ===
 * "viewing"` (or `"viewing-stale"`).
 */
export function InlineEditableText({
  value, onCommit, placeholder, className, inputClassName, testId,
  readOnly = false,
}: Props): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => { if (!editing) setDraft(value); }, [value, editing]);
  useEffect(() => {
    if (editing && inputRef.current) {
      inputRef.current.focus();
      inputRef.current.select();
    }
  }, [editing]);

  if (readOnly) {
    return (
      <span data-testid={testId} className={className}>
        {value === "" ? <span className="text-fg-dim">{placeholder ?? ""}</span> : value}
      </span>
    );
  }

  if (editing) {
    return (
      <span data-interactable="edit">
        <input
          ref={inputRef}
          type="text"
          data-testid={testId}
          value={draft}
          placeholder={placeholder}
          onChange={(e) => setDraft(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter") { onCommit(draft); setEditing(false); }
            else if (e.key === "Escape") { setDraft(value); setEditing(false); }
          }}
          onBlur={() => { onCommit(draft); setEditing(false); }}
          className={inputClassName ?? "bg-transparent border-b border-border outline-none"}
        />
      </span>
    );
  }

  return (
    <span
      data-interactable="edit"
      data-testid={testId}
      onClick={() => setEditing(true)}
      role="button"
      tabIndex={0}
      onKeyDown={(e) => { if (e.key === "Enter") setEditing(true); }}
      className={`${className ?? ""} cursor-text hover:underline decoration-dotted underline-offset-4`}
    >
      {value === "" ? <span className="text-fg-dim">{placeholder ?? ""}</span> : value}
    </span>
  );
}
