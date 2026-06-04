import { useEffect, useState } from "react";
import { Kicker } from "../ui/Kicker";

export interface NotesMarginProps {
  notes: string | null;
  /** Fired on blur with the current draft. */
  onSaveNotes: (notes: string) => void;
  /** Hovering a `q ≈ N.NNN` ref-span emits its numeric q; leaving emits undefined. */
  onHoverQ?: (q: number | undefined) => void;
  /** PLACEMENT ONLY — appended to the aside's className. */
  className?: string;
}

// Matches `q ≈ 0.064` / `q≈0.128` etc. (whitespace around ≈ optional).
const Q_REF = /q\s*≈\s*\d+(?:\.\d+)?/g;
// Capture group to extract the numeric portion.
const Q_CAPTURE = /q\s*≈\s*(\d+(?:\.\d+)?)/;

/**
 * Render a note body into plain + mono `q ≈ N.NNN` reference runs.
 * When `onHoverQ` is provided, each q-ref span wires mouseEnter/mouseLeave
 * to emit the parsed q value or undefined.
 */
function renderNoteBody(
  text: string,
  onHoverQ?: (q: number | undefined) => void,
): JSX.Element[] {
  const out: JSX.Element[] = [];
  let last = 0;
  let key = 0;
  let m: RegExpExecArray | null;
  Q_REF.lastIndex = 0;
  while ((m = Q_REF.exec(text)) !== null) {
    if (m.index > last) out.push(<span key={key++}>{text.slice(last, m.index)}</span>);
    const token = m[0];
    const captureMatch = Q_CAPTURE.exec(token);
    const parsedQ = captureMatch ? Number(captureMatch[1]) : NaN;
    out.push(
      <span
        key={key++}
        data-testid="focus-notes-qref"
        className="font-mono font-semibold text-ink-soft"
        {...(onHoverQ !== undefined
          ? {
              onMouseEnter: () => onHoverQ(parsedQ),
              onMouseLeave: () => onHoverQ(undefined),
            }
          : {})}
      >
        {token}
      </span>,
    );
    last = m.index + token.length;
  }
  if (last < text.length) out.push(<span key={key++}>{text.slice(last)}</span>);
  return out;
}

/**
 * NotesMargin — greenfield port of FocusNotesMargin.
 *
 * Accepts notes as a plain string (+ save callback) rather than a full Sample,
 * so the Focus page can source notes from the correct cache and own persistence.
 * Draft + focus-gate state is internal.
 *
 * The optional `onHoverQ` prop wires each `q ≈ N.NNN` reference span to light
 * the matching peak in the trace (the q-link).
 */
export function NotesMargin({
  notes,
  onSaveNotes,
  onHoverQ,
  className,
}: NotesMarginProps): JSX.Element {
  const [draft, setDraft] = useState(notes ?? "");
  const [focused, setFocused] = useState(false);

  // Focus-gated sync: while the textarea is focused a background refetch or
  // another user's save must NOT clobber a mid-edit draft. Switching samples
  // changes `notes`, which re-syncs while unfocused (focus leaves first).
  useEffect(() => {
    if (!focused) setDraft(notes ?? "");
  }, [notes, focused]);

  const saved = (notes ?? "").trim();
  const hasNote = saved.length > 0;

  return (
    <aside
      data-testid="focus-notes-margin"
      className={
        "flex flex-col gap-[15px] overflow-y-auto border-l border-hair bg-paper px-[19px] py-[22px]" +
        (className !== undefined ? " " + className : "")
      }
    >
      <div className="flex items-baseline gap-[7px]">
        <Kicker as="span" tone="faint">Notes</Kicker>
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
          {renderNoteBody(saved, onHoverQ)}
        </div>
      )}

      {/* `.nm-add` — dashed "✎ Add a note…" affordance wrapping the focus-gated editor. */}
      <div
        data-testid="focus-notes-add"
        className="mt-auto border-b border-dashed border-hair-strong pb-[7px] pt-[6px]
                   text-xs text-ink-faint before:text-accent before:content-['✎_']"
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
