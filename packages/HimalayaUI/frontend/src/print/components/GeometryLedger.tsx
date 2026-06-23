import type { JSX } from "react";
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
}

export interface GeometryLedgerProps {
  rows: GeometryRow[];
  onOverride: (key: string) => void;
  onRevert: (key: string) => void;
  onUndo: () => void;
  canUndo: boolean;
  onRedo: () => void;
  canRedo: boolean;
  discrepancyCount?: number;
  className?: string;
  /** Key of the row currently being inline-edited (undefined = not editing). */
  editingKey?: string;
  /** Controlled draft value for the row being edited. */
  editDraft?: string;
  /** Called when the user types into the inline input. */
  onEditDraftChange?: (value: string) => void;
  /** Called on Enter or blur -- parent decides whether to PATCH. */
  onEditCommit?: () => void;
  /** Called on Escape -- parent cancels the edit, no PATCH. */
  onEditCancel?: () => void;
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
          const isEditing = p.editingKey === r.key;
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
                  value={p.editDraft ?? ""}
                  onValueChange={p.onEditDraftChange ?? (() => {})}
                  aria-label={`Override ${r.label}`}
                  className="w-32"
                  autoFocus
                  onKeyDown={(e) => {
                    if (e.key === "Enter") { e.preventDefault(); p.onEditCommit?.(); }
                    if (e.key === "Escape") { e.preventDefault(); p.onEditCancel?.(); }
                  }}
                  onBlur={() => p.onEditCommit?.()}
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
                    onClick={() => p.onOverride(r.key)}
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
