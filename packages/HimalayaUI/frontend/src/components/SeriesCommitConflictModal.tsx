/**
 * SeriesCommitConflictModal — the series half of the conflict resolution
 * surface (I3.5b). Reuses `ConflictModalShell` (the extracted chrome) and
 * supplies SERIES semantics: a plate-commit 409 carries the server's current
 * `Series` (`isFullSeries`). It reads the SAME `pendingConflict` slot as the
 * comparison `ConflictModal`; each renders only for its own entity kind.
 *
 * COMMIT-ONLY: recipe-save (`PATCH /api/series/:id`) never 409s, so this modal
 * resolves only the plate-commit conflict.
 *
 * Actions:
 *   - Discard mine  — clear the slot + draft, navigate to /series/:id (the
 *                     server's current state).
 *   - Overwrite     — re-flush `useCommitSeriesPlate` with
 *                     `expected_content_hash = current_hash` (the server's new
 *                     hash) over the server's CURRENT plate members.
 *
 * No "Fork" action (forking a series is not part of I3.5b).
 *
 * Double-click guard: a synchronous in-flight ref, mirroring ConflictModal —
 * `overwriteBusy` (= save.isPending) flips async, so two fast clicks could both
 * pass; the ref bails the second before mutate runs.
 */
import { useCallback, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { ConflictModalShell } from "./ConflictModalShell";
import { useCommitSeriesPlate } from "../queries";
import { buildSeriesCommitBody } from "../lib/series/buildSeriesCommitBody";
import { isFullSeries } from "../api";
import type { Series } from "../api";

export function SeriesCommitConflictModal(): JSX.Element | null {
  const conflict           = useAppState((s) => s.pendingConflict);
  const seriesDraft        = useAppState((s) => s.seriesDraft);
  const setPendingConflict = useAppState((s) => s.setPendingConflict);
  const discardSeriesDraft = useAppState((s) => s.discardSeriesDraft);
  const navigate           = useNavigate();
  const commit             = useCommitSeriesPlate();

  // Render ONLY for a series-shaped conflict (samples + state arrays present).
  const rawState = conflict?.current_state ?? null;
  const serverState: Series | null =
    rawState !== null && isFullSeries(rawState) ? rawState : null;
  const serverHash: string | null = conflict?.current_hash ?? null;
  const isOpen = serverState !== null;

  const closeWithoutAction = useCallback(() => {
    setPendingConflict(null);
  }, [setPendingConflict]);

  const handleDiscard = useCallback(() => {
    if (serverState === null) return;
    discardSeriesDraft();
    setPendingConflict(null);
    navigate(`/series/${serverState.id}`);
  }, [serverState, discardSeriesDraft, setPendingConflict, navigate]);

  // Overwrite — re-commit the server's current plate against the server's new
  // hash. We rebuild from `serverState.members` (the canonical plate) so a
  // stale local snapshot can't ride the overwrite. baseHash is overridden with
  // the server's current_hash via a synthetic draft so the commit body carries
  // `expected_content_hash = current_hash`.
  const overwriteInFlightRef = useRef(false);
  const handleOverwrite = useCallback(() => {
    if (serverState === null || serverHash === null || seriesDraft === null) return;
    if (commit.isPending || overwriteInFlightRef.current) return;
    overwriteInFlightRef.current = true;
    const body = buildSeriesCommitBody(
      { ...seriesDraft, baseHash: serverHash },
      serverState.members,
    );
    commit.mutate({ id: serverState.id, ...body });
  }, [serverState, serverHash, seriesDraft, commit]);

  // Overwrite success → clear slot + draft, navigate to the committed series.
  useEffect(() => {
    if (!commit.isSuccess || !overwriteInFlightRef.current) return;
    overwriteInFlightRef.current = false;
    const resp = commit.data as Series | undefined;
    discardSeriesDraft();
    setPendingConflict(null);
    if (resp && typeof resp.id === "number") navigate(`/series/${resp.id}`);
  }, [commit.isSuccess, commit.data, discardSeriesDraft, setPendingConflict, navigate]);

  // Release the guard on a terminal error (e.g. a SECOND 409) so the user can
  // retry; the bridge re-populates `pendingConflict` with the new server state.
  useEffect(() => {
    if (commit.error && overwriteInFlightRef.current) {
      overwriteInFlightRef.current = false;
    }
  }, [commit.error]);

  if (!isOpen || serverState === null) return null;

  return (
    <ConflictModalShell
      open
      heading="Series changed while you were editing"
      subtitle="Someone else committed this series since you entered edit mode. Choose how to resolve."
      serverPanel={{
        label: "Server (current)",
        testId: "conflict-panel-server",
        title: serverState.title,
        memberCount: serverState.members.length,
        description: serverState.description,
        updatedAt: serverState.updated_at,
      }}
      localPanel={{
        label: "Your draft",
        testId: "conflict-panel-local",
        title: seriesDraft?.title ?? "(no title)",
        memberCount: seriesDraft?.recipe.length ?? 0,
        description: seriesDraft?.description ?? null,
        updatedAt: null,
      }}
      onClose={closeWithoutAction}
      onDiscard={handleDiscard}
      discardLabel="Discard my changes"
      onOverwrite={handleOverwrite}
      overwriteBusy={commit.isPending}
    />
  );
}
