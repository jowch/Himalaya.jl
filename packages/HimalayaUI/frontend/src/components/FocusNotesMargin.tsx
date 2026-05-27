import { useEffect, useState } from "react";
import type { Sample } from "../api";

interface Props {
  sample: Sample;
  /** Called on blur with the edited notes value. */
  onSaveNotes: (notes: string) => void;
}

/**
 * FocusNotesMargin — the Notes margin on the focus workspace (the mockup's
 * `.notes-margin`: a quiet ruled column, deliberately not a card).
 *
 * Focus-gated edit (mirrors SampleMetadataCard / QNumInput): the external
 * `sample.notes` is synced into the draft only while the textarea is NOT
 * focused, so a background refetch or another user's save can't clobber a
 * mid-edit draft. Switching samples still works because focus moves out
 * before the new sample renders.
 *
 * The `sample` prop MUST come from the experiment-scoped `useSamples` cache
 * (the one `updateSampleMutator` patches) — see FocusWorkspaceLayout — so the
 * on-blur re-sync reflects the just-saved edit, not a stale corpus value.
 */
export function FocusNotesMargin({ sample, onSaveNotes }: Props): JSX.Element {
  const [draft, setDraft] = useState(sample.notes ?? "");
  const [focused, setFocused] = useState(false);

  useEffect(() => {
    if (!focused) setDraft(sample.notes ?? "");
  }, [sample.id, sample.notes, focused]);

  return (
    <aside data-testid="focus-notes-margin"
           className="flex flex-col gap-2 border-l border-hair bg-paper px-4 py-5">
      <div className="text-meta uppercase tracking-wider text-ink-faint">Notes</div>
      <textarea
        data-testid="focus-notes-input"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onFocus={() => setFocused(true)}
        onBlur={() => {
          setFocused(false);
          onSaveNotes(draft);
        }}
        placeholder="Add a note…"
        className="min-h-[120px] w-full resize-none bg-transparent text-sm
                   text-ink-soft outline-none placeholder:text-ink-faint"
      />
    </aside>
  );
}
