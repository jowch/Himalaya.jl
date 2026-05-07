/**
 * ConflictModal — Compare page 409-content_hash conflict resolution
 * (Plan §Phase 12).
 *
 * Mounted at App.tsx so a save-in-flight survives navigation away from
 * edit mode. Driven by Zustand's `pendingConflict` slot, which
 * `useSaveComparison` populates on `ConflictError`.
 *
 * Three actions:
 *   1. Discard mine     — drop the local draft, navigate to the server's
 *                         current state in review mode, clear the slot.
 *   2. Overwrite mine   — re-submit `saveComparison` with
 *                         `expected_content_hash = current_hash` from the
 *                         conflict body. Each click mints a fresh
 *                         `client_op_id` via `useQueueMutation`. Success
 *                         clears the modal + navigates; a SECOND 409
 *                         re-opens the modal with the newer `current_state`
 *                         (handled implicitly because `useSaveComparison`'s
 *                         onError calls `setPendingConflict` again, and the
 *                         modal reads from the live Zustand slot).
 *   3. Fork to a new    — `startForkDraft` against the SERVER's current
 *      comparison         state (so the fork starts from canonical truth,
 *                         matching `EditOrForkButton`'s semantics from
 *                         Phase 11). Navigates to /experiments/:eid/compare/new.
 *
 * Esc / outside click: close the modal (clear the slot) WITHOUT committing
 * any of the three actions. The local draft is preserved — Discard is the
 * explicit "I give up" path; closing without acting just dismisses the
 * conflict prompt and returns the user to their edit session.
 *
 * The :eid for navigation is parsed from the current URL pathname. The
 * conflict can only be raised from a save initiated on an experiment-
 * scoped edit route (/experiments/:eid/compare/...) or the experiment-
 * less /compare/all flow; the eid resolver handles both.
 */
import { useCallback, useEffect, useMemo, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { useFocusTrap } from "../hooks/useFocusTrap";
import { useSaveComparison } from "../queries";
import { computeMemberSnapshot } from "../lib/comparison/snapshot";
import type {
  Comparison, ComparisonMemberInput, SaveComparisonBody,
} from "../api";
import type { ActiveDraft } from "../lib/comparison/draft";

/** Extract `:eid` from a `/experiments/:eid/...` URL or return undefined. */
function extractEid(pathname: string): number | undefined {
  const m = pathname.match(/^\/experiments\/(\d+)\b/);
  if (!m) return undefined;
  const n = Number(m[1]);
  return Number.isFinite(n) ? n : undefined;
}

/**
 * Build the SaveComparisonBody from a draft, parameterized by the
 * `expected_content_hash` to use. Used by the "Overwrite" action to
 * re-submit the local draft against the server's NEW current_hash.
 *
 * Mirrors `ComparePageEdit::handleSave` snapshot-recompute semantics so a
 * stale snapshot doesn't ride the overwrite request.
 */
function buildOverwritePayload(
  draft: ActiveDraft,
  serverHash: string,
  qc: ReturnType<typeof useQueryClient>,
): SaveComparisonBody & { id?: number } {
  const members: ComparisonMemberInput[] = draft.members.map((m) => {
    const snapshot = m.exposure_id !== null
      ? computeMemberSnapshot(m.exposure_id, qc)
      : (m.snapshot ?? {
          effective_peaks: [],
          confirmed_index: null,
          analysis_inputs_hash: "",
        });
    const out: ComparisonMemberInput = {
      exposure_id: m.exposure_id,
      display_order: m.display_order,
      band_height: m.band_height,
      y_offset: m.y_offset,
      normalization: m.normalization,
      snapshot,
    };
    if (m.id !== undefined) out.id = m.id;
    if (m.color_override !== undefined) out.color_override = m.color_override;
    if (m.label_override !== undefined) out.label_override = m.label_override;
    if (m.q_window_min !== undefined) out.q_window_min = m.q_window_min;
    if (m.q_window_max !== undefined) out.q_window_max = m.q_window_max;
    if (m.peak_display !== undefined) out.peak_display = m.peak_display;
    return out;
  });
  const payload: SaveComparisonBody & { id?: number } = {
    title: draft.title,
    members,
    expected_content_hash: serverHash,
  };
  if (draft.id !== undefined) payload.id = draft.id;
  if (draft.description !== "") payload.description = draft.description;
  if (draft.forkedFromId !== undefined) payload.forked_from_id = draft.forkedFromId;
  if (draft.forkedAtHash !== undefined) payload.forked_at_hash = draft.forkedAtHash;
  return payload;
}

export function ConflictModal(): JSX.Element | null {
  const conflict           = useAppState((s) => s.pendingConflict);
  const draft              = useAppState((s) => s.activeDraft);
  const setPendingConflict = useAppState((s) => s.setPendingConflict);
  const discardDraft       = useAppState((s) => s.discardDraft);
  const startForkDraft     = useAppState((s) => s.startForkDraft);
  const navigate           = useNavigate();
  const location           = useLocation();
  const qc                 = useQueryClient();
  const save               = useSaveComparison();

  const dialogRef = useRef<HTMLDivElement>(null);
  const isOpen    = conflict !== null;
  useFocusTrap(dialogRef, isOpen);

  const serverState: Comparison | null = conflict?.current_state ?? null;
  const serverHash:  string | null     = conflict?.current_hash  ?? null;

  const eid = useMemo(() => extractEid(location.pathname), [location.pathname]);

  const goToReview = useCallback(
    (comparisonId: number) => {
      if (eid !== undefined) navigate(`/experiments/${eid}/compare/${comparisonId}`);
      else                   navigate("/compare/all");
    },
    [navigate, eid],
  );

  const goToNewDraft = useCallback(() => {
    if (eid !== undefined) navigate(`/experiments/${eid}/compare/new`);
    else                   navigate("/compare/all");
  }, [navigate, eid]);

  // Esc and outside-click both close the modal WITHOUT committing — the
  // local draft is preserved. This matches the design call (non-destructive
  // dismiss); Discard is the explicit "I give up" path.
  const closeWithoutAction = useCallback(() => {
    setPendingConflict(null);
  }, [setPendingConflict]);

  // Esc handler — bound at document level so it works even when focus is
  // briefly outside the dialog (e.g. after a focus restore between renders).
  useEffect(() => {
    if (!isOpen) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") {
        e.preventDefault();
        closeWithoutAction();
      }
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [isOpen, closeWithoutAction]);

  // Discard: drop our draft, go to the server's current state, clear slot.
  const handleDiscard = useCallback(() => {
    if (serverState === null) return;
    discardDraft();
    setPendingConflict(null);
    goToReview(serverState.id);
  }, [serverState, discardDraft, setPendingConflict, goToReview]);

  // Overwrite: re-submit with expected_content_hash = server's current_hash.
  // The `useSaveComparison` post-success effect lives in ComparePageEdit, but
  // the modal navigates itself when the mutation succeeds — `save.isSuccess`
  // flips to true and we tear down the slot + navigate.
  const overwriteInFlightRef = useRef(false);
  const handleOverwrite = useCallback(() => {
    if (draft === null || serverHash === null) return;
    overwriteInFlightRef.current = true;
    save.mutate(buildOverwritePayload(draft, serverHash, qc));
  }, [draft, serverHash, save, qc]);

  // Watch for overwrite success → clear modal + navigate.
  // We deliberately read `save.data` rather than the conflict slot because the
  // save mutator's success bumps `data` to the canonical Comparison response.
  // A SECOND 409 leaves `isSuccess` false and (re-)sets `pendingConflict` via
  // `useSaveComparison`'s effect — the modal stays open with new server state.
  useEffect(() => {
    if (!save.isSuccess) return;
    if (!overwriteInFlightRef.current) return;
    overwriteInFlightRef.current = false;
    const response = save.data;
    if (response && typeof response.id === "number") {
      discardDraft();
      setPendingConflict(null);
      goToReview(response.id);
    }
  }, [save.isSuccess, save.data, discardDraft, setPendingConflict, goToReview]);

  // Fork: start a fork-flavored draft from the SERVER's current state (so
  // the fork inherits the canonical truth, not the user's stale-baseHash
  // local edits). Mirrors `EditOrForkButton::onFork` from Phase 11.
  // Documenting the shape choice: forking from the server's current state
  // matches the Phase 11 EditOrForkButton; forking from the local draft
  // would silently lose the conflict context (the user's local title edit
  // would survive but the lineage would still point at the server hash).
  // Server-state-as-fork-source is the least surprising: the user explicitly
  // chose "Fork", which traditionally means "branch off the canonical state".
  const handleFork = useCallback(() => {
    if (serverState === null) return;
    startForkDraft(serverState, qc);
    setPendingConflict(null);
    goToNewDraft();
  }, [serverState, startForkDraft, qc, setPendingConflict, goToNewDraft]);

  if (!isOpen || serverState === null) return null;

  const localMemberCount  = draft?.members.length ?? 0;
  const serverMemberCount = serverState.members.length;
  const localTitle  = draft?.title ?? "(no title)";
  const serverTitle = serverState.title;

  return (
    <div
      data-testid="conflict-modal"
      className="fixed inset-0 z-50 flex items-center justify-center
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm
                 anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) closeWithoutAction(); }}
    >
      <div
        ref={dialogRef}
        role="dialog"
        aria-modal="true"
        aria-labelledby="conflict-title"
        aria-describedby="conflict-subtitle"
        className="w-[min(820px,calc(100vw-48px))] max-h-[80vh]
                   bg-bg-elevated border border-border rounded-xl shadow-2xl
                   flex flex-col overflow-hidden anim-pal-scale"
      >
        <header className="px-5 py-4 border-b border-border">
          <h2 id="conflict-title" className="text-fg text-lg font-medium">
            Comparison changed while you were editing
          </h2>
          <p id="conflict-subtitle" className="text-fg-muted text-sm mt-1">
            Someone else submitted updates since you entered edit mode.
            Choose how to resolve.
          </p>
        </header>

        <div className="flex-1 min-h-0 overflow-y-auto grid grid-cols-2 gap-3 p-5">
          <Panel
            label="Server (current)"
            testId="conflict-panel-server"
            title={serverTitle}
            memberCount={serverMemberCount}
            description={serverState.description}
            updatedAt={serverState.updated_at}
          />
          <Panel
            label="Your draft"
            testId="conflict-panel-local"
            title={localTitle}
            memberCount={localMemberCount}
            description={draft?.description ?? null}
            updatedAt={null}
          />
        </div>

        <footer className="flex items-center gap-2 px-5 py-3 border-t border-border">
          <button
            type="button"
            data-testid="conflict-discard"
            onClick={handleDiscard}
            className="px-3 py-1.5 rounded border border-border text-fg text-sm
                       hover:bg-bg-hover"
          >
            Discard my changes
          </button>
          <button
            type="button"
            data-testid="conflict-fork"
            onClick={handleFork}
            className="px-3 py-1.5 rounded border border-border text-fg text-sm
                       hover:bg-bg-hover"
          >
            Fork to a new comparison
          </button>
          <span className="flex-1" />
          <button
            type="button"
            data-testid="conflict-overwrite"
            onClick={handleOverwrite}
            disabled={save.isPending}
            className="px-3 py-1.5 rounded border border-accent bg-accent
                       text-bg text-sm disabled:opacity-60"
          >
            {save.isPending ? "Saving…" : "Overwrite with mine"}
          </button>
        </footer>
      </div>
    </div>
  );
}

interface PanelProps {
  label: string;
  testId: string;
  title: string;
  memberCount: number;
  description: string | null;
  updatedAt: string | null;
}

function Panel({
  label, testId, title, memberCount, description, updatedAt,
}: PanelProps): JSX.Element {
  return (
    <section
      data-testid={testId}
      className="border border-border rounded-md p-3 flex flex-col gap-2 min-w-0"
    >
      <header className="text-xs uppercase tracking-wide text-fg-dim">
        {label}
      </header>
      <div className="text-fg font-medium truncate" data-testid={`${testId}-title`}>
        {title || "(no title)"}
      </div>
      <div className="text-fg-muted text-sm" data-testid={`${testId}-members`}>
        {memberCount} {memberCount === 1 ? "member" : "members"}
      </div>
      {description && (
        <div
          className="text-fg-muted text-sm whitespace-pre-wrap"
          data-testid={`${testId}-description`}
        >
          {description}
        </div>
      )}
      {updatedAt && (
        <div className="text-fg-dim text-xs mt-auto" data-testid={`${testId}-updated`}>
          Updated {updatedAt}
        </div>
      )}
    </section>
  );
}
