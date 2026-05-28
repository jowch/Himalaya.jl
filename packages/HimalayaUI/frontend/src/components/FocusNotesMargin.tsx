import { useEffect, useState } from "react";
import type { Sample } from "../api";

interface Props {
  sample: Sample;
  /** Called on blur with the edited notes value. */
  onSaveNotes: (notes: string) => void;
}

// Matches the mockup's mono `q ≈ 0.064` reference (whitespace around ≈ optional,
// so a "q≈0.064" note also formats). Stage-1 detection only — no schema change.
const Q_REF = /q\s*≈\s*\d+(?:\.\d+)?/g;

/**
 * Split a note body into plain + mono `q ≈ N.NNN` reference runs so the
 * references mirror the q-link triple's monospace voice (DESIGN.md §3
 * Monospace-Means-Measurement). Returns React nodes, not HTML.
 */
function renderNoteBody(text: string): JSX.Element[] {
  const out: JSX.Element[] = [];
  let last = 0;
  let key = 0;
  let m: RegExpExecArray | null;
  Q_REF.lastIndex = 0;
  while ((m = Q_REF.exec(text)) !== null) {
    if (m.index > last) out.push(<span key={key++}>{text.slice(last, m.index)}</span>);
    out.push(
      <span
        key={key++}
        data-testid="focus-notes-qref"
        className="font-mono font-semibold text-ink-soft"
      >
        {m[0]}
      </span>,
    );
    last = m.index + m[0].length;
  }
  if (last < text.length) out.push(<span key={key++}>{text.slice(last)}</span>);
  return out;
}

/**
 * FocusNotesMargin — the Notes margin on the focus workspace (the mockup's
 * `.notes-margin`: a quiet ruled column, deliberately not a card).
 *
 * Stage-1 Print treatment (R3-F04): a `Notes` head + count badge, the saved
 * note rendered as quiet body text with mono-formatted `q ≈ N.NNN` references,
 * and a dashed `✎ Add a note…` affordance wrapping the editor. The bare
 * textarea-first layout is gone. Threaded per-note metadata (author initials +
 * timestamp, mockup `.note .n-meta`) is stage 2, deferred until the schema
 * carries per-note rows — this stage renders the single `sample.notes` field.
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

  const saved = (sample.notes ?? "").trim();
  const hasNote = saved.length > 0;

  return (
    <aside
      data-testid="focus-notes-margin"
      className="flex flex-col gap-[15px] overflow-y-auto border-l border-hair bg-paper px-[19px] py-[22px]"
    >
      <div className="flex items-baseline gap-[7px]">
        <span className="text-meta uppercase tracking-wider text-ink-faint">Notes</span>
        {hasNote && (
          <span
            data-testid="focus-notes-count"
            className="font-mono text-meta font-bold text-ink-faint"
          >
            1
          </span>
        )}
      </div>

      {hasNote && (
        <div
          data-testid="focus-notes-body"
          className="text-xs leading-[1.6] text-ink-soft"
        >
          {renderNoteBody(saved)}
        </div>
      )}

      {/* `.nm-add` — the dashed "✎ Add a note…" affordance from the mockup,
          wrapping the focus-gated editor so the textarea no longer reads as a
          bare academic-utility field. */}
      <div
        data-testid="focus-notes-add"
        className="mt-auto border-b border-dashed border-hair-strong pb-[7px] pt-[6px]
                   text-xs text-ink-faint before:text-[var(--color-accent)] before:content-['✎_']"
      >
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
          className="min-h-[44px] w-full resize-none bg-transparent text-xs
                     text-ink-soft outline-none placeholder:text-ink-faint"
        />
      </div>
    </aside>
  );
}
