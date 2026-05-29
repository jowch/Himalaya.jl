import { useState } from "react";

interface Props {
  sampleId: number;
  sampleName: string;
  value: string;
  flagged: boolean;
  onChange: (value: string) => void; // commit a new value (also clears the flag)
  onToggleFlag: () => void; // accept the read / re-flag (undo-tracked by parent)
}

/**
 * Confirm-not-fill-out value cell (series-scoping.html `.s-val`, finding S-E).
 *  - Confident: ink text; click re-opens it as an input (every value stays
 *    re-openable, including one resolved by mistake).
 *  - Flagged: amber dashed underline + "▸ check the read"; click accepts the
 *    read (clears the flag) — never a permanent text input.
 *  - Editing: a transient input; blur/Enter commits via onChange, Escape cancels.
 */
export function ScopingValueCell({
  sampleId, sampleName, value, flagged, onChange, onToggleFlag,
}: Props): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);

  function open(): void {
    if (flagged) { onToggleFlag(); return; } // accept the read in one click
    setDraft(value);
    setEditing(true);
  }
  function commit(): void {
    setEditing(false);
    if (draft !== value) onChange(draft);
  }

  if (editing) {
    return (
      <input
        type="text"
        data-testid={`scoping-value-input-${sampleId}`}
        aria-label={`Value for ${sampleName}`}
        value={draft}
        autoFocus
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => {
          if (e.key === "Enter") commit();
          if (e.key === "Escape") { setDraft(value); setEditing(false); }
        }}
        className="w-24 shrink-0 rounded border border-hair-strong bg-plate px-2 py-1 text-right font-mono text-base text-ink"
      />
    );
  }

  return (
    <button
      type="button"
      data-testid={`scoping-value-${sampleId}`}
      data-flagged={flagged ? "true" : undefined}
      onClick={open}
      title={flagged ? "Click to accept this read" : "Click to re-open this value"}
      className={`group shrink-0 text-right ${flagged ? "text-print-accent" : "text-ink"}`}
    >
      <span
        className={`font-mono text-base font-bold ${
          flagged
            ? "border-b-[1.5px] border-dashed border-print-accent/60 pb-px"
            : "border-b-[1.5px] border-dotted border-transparent pb-px group-hover:border-hair-strong"
        }`}
      >
        {value || "—"}
      </span>
      {flagged ? (
        <span className="mt-0.5 block text-xs font-bold uppercase tracking-wide text-print-accent">
          ▸ check the read
        </span>
      ) : null}
    </button>
  );
}
