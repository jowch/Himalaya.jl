import type { JSX } from "react";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";

export type GeometrySource = "prp" | "setup" | "user" | "default";

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
  discrepancyCount?: number;
  className?: string;
}

// Source label display strings -- no em dashes.
const SOURCE_LABEL: Record<GeometrySource, string> = {
  prp: "PRP",
  setup: "setup files",
  user: "edited",
  default: "unset",
};

// SourceChip: a small presentational badge for the provenance of a geometry
// field. "edited" (user) gets the accent-wash tint; all others get a neutral
// paper-sunk background. This component lives in print/components/ (NOT the
// design-guard-exempt print/ui/ layer) so it may only use named token utilities
// -- no bg-[...] or color-mix literals. bg-accent-wash is sanctioned because
// --color-accent-wash is declared in styles.css @theme.
function SourceChip({ source }: { source: GeometrySource }): JSX.Element {
  const isUser = source === "user";
  return (
    <span
      className={
        isUser
          ? "rounded-sm bg-accent-wash px-1.5 py-0.5 text-xs font-bold uppercase text-accent"
          : "rounded-sm bg-paper-sunk px-1.5 py-0.5 text-xs font-bold uppercase text-ink-faint"
      }
    >
      {SOURCE_LABEL[source]}
    </span>
  );
}

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
        {p.rows.map((r) => (
          <div
            key={r.key}
            className="flex items-center gap-3 border-b border-hair py-2.5 last:border-b-0"
          >
            <span className="w-32 shrink-0 text-xs text-ink-soft">{r.label}</span>
            <span className="font-mono text-sm font-medium text-ink">{r.value}</span>
            <SourceChip source={r.source} />
            {r.source === "user" ? (
              <Button
                variant="ghost"
                aria-label={`Revert ${r.label}`}
                onClick={() => p.onRevert(r.key)}
              >
                Revert
              </Button>
            ) : null}
            <Button
              variant="ghost"
              aria-label={`Override ${r.label}`}
              onClick={() => p.onOverride(r.key)}
            >
              {r.source === "user" ? "Edit" : "Override"}
            </Button>
          </div>
        ))}
      </div>
    </Card>
  );
}
