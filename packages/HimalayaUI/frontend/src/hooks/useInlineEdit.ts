import { useEffect, useRef, useState, type RefObject } from "react";

export interface InlineEdit<K> {
  /** The row/field currently being edited, or null. Boolean callers pass a
   *  sentinel key (e.g. `"_"`). */
  editingKey: K | null;
  /** Current draft text. */
  draft: string;
  setDraft: (v: string) => void;
  /** Attach to the edit `<input>`; focused + selected when editing opens. */
  inputRef: RefObject<HTMLInputElement>;
  /** Enter edit mode for `key`, seeding the draft (and the changed-check) with `seed`. */
  begin: (key: K, seed: string) => void;
  /** Commit: trims, fires `onCommit` only if the value actually changed, then
   *  closes. Guarded so a blur firing right after Enter cannot commit twice. */
  commit: () => void;
  /** Leave edit mode without committing. */
  cancel: () => void;
}

/**
 * One inline edit-in-place state machine for every "click a value, type, Enter"
 * field (sample rename, experiment description, scoping value, file patterns,
 * geometry params). Unified semantics: focus+select on open, trim on commit,
 * skip no-op commits, and a double-commit guard.
 *
 * Field-specific behaviour stays in the caller, NOT in this hook:
 *  - what an empty value means (ignore vs clear-to-null) → handle in `onCommit`;
 *  - the commit key (Enter vs Cmd+Enter for multiline) → wire in the caller's
 *    own `onKeyDown`, calling `commit()`.
 */
export function useInlineEdit<K>(
  onCommit: (key: K, value: string) => void,
): InlineEdit<K> {
  const [editingKey, setEditingKey] = useState<K | null>(null);
  const [draft, setDraft] = useState("");
  const seedRef = useRef("");
  // Synchronous guard: a blur firing right after Enter must not commit twice (a
  // double commit would push two undo entries for one edit). State alone is too
  // late — the blur handler closes over the pre-update render — so the ref is
  // cleared inside commit() before the second call can read it.
  const guardRef = useRef<K | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (editingKey !== null) {
      inputRef.current?.focus();
      inputRef.current?.select();
    }
  }, [editingKey]);

  const begin = (key: K, seed: string): void => {
    seedRef.current = seed;
    guardRef.current = key;
    setEditingKey(key);
    setDraft(seed);
  };

  const commit = (): void => {
    const key = guardRef.current;
    if (key === null) return;
    guardRef.current = null;
    setEditingKey(null);
    const next = draft.trim();
    if (next !== seedRef.current) onCommit(key, next);
  };

  const cancel = (): void => {
    guardRef.current = null;
    setEditingKey(null);
  };

  return { editingKey, draft, setDraft, inputRef, begin, commit, cancel };
}
