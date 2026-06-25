import { type JSX } from "react";
import { useInlineEdit } from "../../hooks/useInlineEdit";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";
import { SourceChip } from "./SourceChip";

export type GeometrySource = "prp" | "setup" | "user" | "default" | "computed";

export interface GeometryRow {
  key: string;
  label: string;
  value: string;
  source: GeometrySource;
  /** Raw, unit-less seed for the inline editor (e.g. "9.00", not "9.00 keV").
   *  Undefined means the row has no editable value yet -- Override is a no-op. */
  editValue?: string | undefined;
}

export interface GeometryLedgerProps {
  rows: GeometryRow[];
  /** Commit an inline override: fired on Enter/blur with the row key and the
   *  raw draft string. The parent parses it, decides no-op vs PATCH, and owns
   *  the undo stack -- this component only owns which row is open and its draft. */
  onCommit: (key: string, value: string) => void;
  onRevert: (key: string) => void;
  onUndo: () => void;
  canUndo: boolean;
  onRedo: () => void;
  canRedo: boolean;
  discrepancyCount?: number;
  className?: string;
}

// Source label display strings -- no em dashes.
const SOURCE_LABEL: Record<GeometrySource, string> = {
  prp: "PRP",
  setup: "setup files",
  user: "edited",
  default: "unset",
  computed: "computed",
};

export function GeometryLedger(p: GeometryLedgerProps): JSX.Element {
  // The hook owns the blur-after-Enter double-commit guard and skips no-ops; the
  // parent still parses the (now trimmed) value and decides no-op vs PATCH.
  const { editingKey, draft, setDraft, inputRef, begin, commit, cancel } =
    useInlineEdit<string>(p.onCommit);

  return (
    <Card className={p.className}>
      <div className="flex items-center justify-between border-b border-hair px-4 py-3.5">
        <h3 className="text-headline text-ink">Geometry</h3>
        <div className="flex items-center gap-3.5">
          {p.canUndo ? (
            <button
              type="button"
              className="text-xs font-semibold text-accent"
              onClick={p.onUndo}
              aria-label="Undo last change"
            >
              Undo last change
            </button>
          ) : null}
          {p.canRedo ? (
            <button
              type="button"
              className="text-xs font-semibold text-accent"
              onClick={p.onRedo}
              aria-label="Redo last change"
            >
              Redo
            </button>
          ) : null}
          <span className="text-xs text-ink-faint">auto-derived</span>
        </div>
      </div>

      {p.discrepancyCount !== undefined && p.discrepancyCount > 0 ? (
        <div className="m-4 rounded-sm border border-warning bg-paper-sunk px-3.5 py-2.5 text-xs text-ink-soft">
          Geometry check found {p.discrepancyCount} issue
          {p.discrepancyCount > 1 ? "s" : ""} -- constant fields varied across PRPs.
        </div>
      ) : null}

      <div className="px-4 pb-3.5">
        {p.rows.map((r) => {
          const isEditing = editingKey === r.key;
          return (
            <div
              key={r.key}
              className="flex items-center gap-3 border-b border-hair py-2.5 last:border-b-0"
            >
              <span className="w-32 shrink-0 text-xs text-ink-soft">{r.label}</span>

              {/* Fixed-width value column so the source chip + actions line up
                  in aligned columns across every row (no ragged float). */}
              {isEditing ? (
                <Input
                  mono
                  inputSize="sm"
                  value={draft}
                  onValueChange={setDraft}
                  aria-label={`Override ${r.label}`}
                  className="w-32"
                  inputRef={inputRef}
                  onKeyDown={(e) => {
                    if (e.key === "Enter") { e.preventDefault(); commit(); }
                    if (e.key === "Escape") { e.preventDefault(); cancel(); }
                  }}
                  onBlur={() => commit()}
                />
              ) : (
                <span className="w-32 shrink-0 font-mono text-sm font-medium text-ink">{r.value}</span>
              )}

              <SourceChip label={SOURCE_LABEL[r.source]} emphasized={r.source === "user"} />

              <div className="ml-auto flex items-center gap-2">
                {r.source === "user" ? (
                  <Button
                    variant="ghost"
                    aria-label={`Revert ${r.label}`}
                    onClick={() => p.onRevert(r.key)}
                  >
                    Revert
                  </Button>
                ) : null}
                {!isEditing && (
                  <Button
                    variant="ghost"
                    aria-label={`Override ${r.label}`}
                    onClick={() => {
                      if (r.editValue !== undefined) begin(r.key, r.editValue);
                    }}
                  >
                    {r.source === "user" ? "Edit" : "Override"}
                  </Button>
                )}
              </div>
            </div>
          );
        })}
      </div>
    </Card>
  );
}
