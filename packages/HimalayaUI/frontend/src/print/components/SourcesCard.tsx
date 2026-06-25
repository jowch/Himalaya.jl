import { useState, type JSX } from "react";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";

export interface SourceRow {
  key: string;
  label: string;
  value: string;
  /** Whether this row's value can be edited in place. Defaults to false.
   *  Directories (data_dir, analysis_dir) are read-only; file patterns are
   *  editable and trigger a rescan when changed. */
  editable?: boolean;
}

export interface SourcesCardProps {
  rows: SourceRow[];
  /** Called when the user commits an edit. Only fires for rows with editable=true. */
  onEdit: (key: string, value: string) => void;
  /** Called when the user clicks "Rescan now". */
  onRescan: () => void;
  className?: string;
}

/** Configuration card showing data directories (read-only) and file patterns
 *  (editable in place). Editing a pattern value and confirming calls onEdit
 *  with the row key and new value. A footer provides the Rescan now action.
 *
 *  Presentational: the page owns the row data, mutation wiring, and onRescan. */
export function SourcesCard(p: SourcesCardProps): JSX.Element {
  const [editing, setEditing] = useState<string | null>(null);
  const [draft, setDraft] = useState("");

  const begin = (r: SourceRow) => {
    setEditing(r.key);
    setDraft(r.value);
  };

  const commit = () => {
    if (editing) {
      p.onEdit(editing, draft.trim());
    }
    setEditing(null);
  };

  const cancel = () => setEditing(null);

  return (
    <Card className={p.className}>
      <div className="flex items-center justify-between border-b border-hair px-4 py-3.5">
        <h3 className="text-headline text-ink">Sources</h3>
        <span className="text-xs text-ink-faint">edits apply on the next rescan</span>
      </div>

      <div className="px-4 pb-2">
        {p.rows.map((r) => (
          <div
            key={r.key}
            className="flex items-center gap-3 border-b border-hair py-2.5 last:border-b-0"
          >
            <span className="w-40 shrink-0 text-xs text-ink-soft">{r.label}</span>

            {editing === r.key ? (
              <Input
                mono
                autoFocus
                value={draft}
                onValueChange={setDraft}
                aria-label={r.label}
                onKeyDown={(e) => {
                  if (e.key === "Enter") commit();
                  if (e.key === "Escape") cancel();
                }}
                onBlur={commit}
              />
            ) : r.editable ? (
              // Editable pattern rows: the value is a button that switches to an
              // input on click. No inline appearance utilities -- text-ink +
              // font-mono are design-token classes, not raw utilities.
              <button
                type="button"
                className="font-mono text-sm text-ink"
                onClick={() => begin(r)}
              >
                {r.value}
              </button>
            ) : (
              // Read-only directory rows: plain span, no click affordance.
              <span className="font-mono text-sm text-ink-soft">{r.value}</span>
            )}
          </div>
        ))}
      </div>

      <div className="flex items-center justify-between gap-4 border-t border-hair px-4 py-3.5">
        <span className="max-w-[64ch] text-xs text-ink-soft">
          Rescan checks for new files and ingests only those. It never re-reads
          existing exposures or overwrites your edits.
        </span>
        <Button variant="outline" onClick={p.onRescan}>
          Rescan now
        </Button>
      </div>
    </Card>
  );
}
