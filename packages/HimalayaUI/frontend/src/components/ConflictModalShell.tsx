/**
 * ConflictModalShell — the presentational chrome shared by the comparison
 * `ConflictModal` and the series `SeriesCommitConflictModal` (I3.5b
 * chrome-extraction). It owns ONLY the dialog frame, focus trap, Esc /
 * outside-click dismiss, the two side-by-side state panels, and the footer
 * action buttons (with the synchronous double-click guard on Overwrite).
 *
 * It knows nothing about comparisons / series / drafts / mutations: the parent
 * supplies the panel contents and the action handlers. This keeps the dense,
 * deeply-coupled semantics out of the shell so reuse can't regress either path.
 *
 * Double-click guard: `disabled={overwriteBusy}` flips ASYNC after the first
 * click, so a fast double-click can pass the disabled check twice. The parent
 * owns the synchronous in-flight ref (each path's overwrite mints a fresh
 * client_op_id, so the idempotency layer can't dedupe a racing pair); the
 * shell just renders `overwriteBusy` and forwards the click.
 */
import { useEffect, useRef } from "react";
import type { ReactNode } from "react";
import { useFocusTrap } from "../hooks/useFocusTrap";

export interface ConflictPanelData {
  label: string;
  testId: string;
  title: string;
  memberCount: number;
  description: string | null;
  updatedAt: string | null;
}

export interface ConflictModalShellProps {
  open: boolean;
  /** Heading + subtitle copy (differs per entity). */
  heading: string;
  subtitle: string;
  serverPanel: ConflictPanelData;
  localPanel: ConflictPanelData;
  /** Esc / outside-click — dismiss without committing (draft preserved). */
  onClose: () => void;
  onDiscard: () => void;
  discardLabel: string;
  onOverwrite: () => void;
  overwriteBusy: boolean;
  /** Optional extra footer action between Discard and the spacer (e.g. Fork). */
  extraAction?: ReactNode;
}

export function ConflictModalShell({
  open, heading, subtitle, serverPanel, localPanel,
  onClose, onDiscard, discardLabel, onOverwrite, overwriteBusy, extraAction,
}: ConflictModalShellProps): JSX.Element | null {
  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, open);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        e.preventDefault();
        onClose();
      }
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [open, onClose]);

  if (!open) return null;

  return (
    <div
      data-testid="conflict-modal"
      className="fixed inset-0 z-50 flex items-center justify-center
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm
                 anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div
        ref={dialogRef}
        role="dialog"
        aria-modal="true"
        aria-labelledby="conflict-title"
        aria-describedby="conflict-subtitle"
        className="w-[min(820px,calc(100vw-48px))] max-h-[80vh]
                   bg-plate border border-hair-strong rounded-xl shadow-2xl
                   flex flex-col overflow-hidden anim-pal-scale"
      >
        <header className="px-5 py-4 border-b border-hair-strong">
          <h2 id="conflict-title" className="text-ink text-lg font-medium">
            {heading}
          </h2>
          <p id="conflict-subtitle" className="text-ink-soft text-sm mt-1">
            {subtitle}
          </p>
        </header>

        <div className="flex-1 min-h-0 overflow-y-auto grid grid-cols-2 gap-3 p-5">
          <Panel {...serverPanel} />
          <Panel {...localPanel} />
        </div>

        <footer className="flex items-center gap-2 px-5 py-3 border-t border-hair-strong">
          <button
            type="button"
            data-testid="conflict-discard"
            onClick={onDiscard}
            className="px-3 py-1.5 rounded border border-hair-strong text-ink text-sm
                       hover:bg-paper-sunk"
          >
            {discardLabel}
          </button>
          {extraAction}
          <span className="flex-1" />
          <button
            type="button"
            data-testid="conflict-overwrite"
            onClick={onOverwrite}
            disabled={overwriteBusy}
            className="px-3 py-1.5 rounded border border-accent bg-accent
                       text-paper text-sm disabled:opacity-60"
          >
            {overwriteBusy ? "Saving…" : "Overwrite with mine"}
          </button>
        </footer>
      </div>
    </div>
  );
}

function Panel({
  label, testId, title, memberCount, description, updatedAt,
}: ConflictPanelData): JSX.Element {
  return (
    <section
      data-testid={testId}
      className="border border-hair-strong rounded-md p-3 flex flex-col gap-2 min-w-0"
    >
      <header className="text-xs uppercase tracking-wide text-ink-faint">
        {label}
      </header>
      <div className="text-ink font-medium truncate" data-testid={`${testId}-title`}>
        {title || "(no title)"}
      </div>
      <div className="text-ink-soft text-sm" data-testid={`${testId}-members`}>
        {memberCount} {memberCount === 1 ? "member" : "members"}
      </div>
      {description && (
        <div
          className="text-ink-soft text-sm whitespace-pre-wrap"
          data-testid={`${testId}-description`}
        >
          {description}
        </div>
      )}
      {updatedAt && (
        <div className="text-ink-faint text-xs mt-auto" data-testid={`${testId}-updated`}>
          Updated {updatedAt}
        </div>
      )}
    </section>
  );
}
