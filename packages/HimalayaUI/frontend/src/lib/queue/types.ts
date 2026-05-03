import type { QueryClient } from "@tanstack/react-query";

// ---------------------------------------------------------------------------
// Optimistic-id invariant
// ---------------------------------------------------------------------------
// Mutators that need an entity id before the server has assigned one
// (e.g. peak_added inserting a peak_curations row, speculative_created
// inserting an indices row) use NEGATIVE placeholder ids in optimistic cache
// writes. Real DB ids are always positive (INTEGER PRIMARY KEY AUTOINCREMENT
// in SQLite never returns ≤ 0 for fresh inserts). Consumers that read these
// caches (PeakRow, MentionChip, etc.) must tolerate negative ids: do not
// `id > 0` filter, do not strict-parse, do not dereference into URLs without
// a sign check. The placeholder is replaced with the real id when the
// mutator's onSuccess runs against the server response.
//
// This is a load-bearing invariant — breaking it surfaces as flicker /
// duplicate rows during optimistic-confirm transitions.

/**
 * The set of mutation kinds the queue framework can carry. Each kind maps to
 * a route handler on the server and (typically) an `apply_event!` event kind.
 */
export type OpKind =
  | "peak_added" | "peak_excluded" | "peak_unexcluded" | "peak_removed"
  | "index_confirmed" | "index_unconfirmed"
  | "speculative_created" | "speculative_deleted"
  | "set_exposure_status" | "select_exposure"
  | "add_tag" | "remove_tag"
  | "post_message" | "update_sample"
  | "reanalyze_exposure"
  | "delete_index";

/**
 * A queued operation: its kind, its per-call client_op_id (Stripe-style
 * idempotency token, fresh per mutate() call), the mutator-specific payload,
 * and an optional exposureId for hook-level scoping.
 */
export interface OpPayload<T = unknown> {
  kind: OpKind;
  payload: T;
  clientOpId: string;
  exposureId?: number;
}

/**
 * The contract returned by a mutator's `onMutate` callback. Holds enough
 * information to revert the optimistic cache write on rollback (replay
 * coordinator) or on HTTP error.
 */
export interface RollbackContext {
  restore: () => void;
}

/**
 * Deferred-promise: standard pattern for resolving a Promise from outside
 * its constructor. Used by the replay coordinator to resolve a pending
 * mutation when its corresponding SSE event arrives.
 */
export interface PendingDeferred<T> {
  promise: Promise<T>;
  resolve: (value: T) => void;
  reject: (reason: unknown) => void;
}

/**
 * The shape of an SSE frame as parsed from the JSON `data:` line. Mirrors
 * the Julia-side `broadcast_event!` JSON shape (events.jl).
 *
 * Optional `post_state` carries an enriched snapshot for replay-without-
 * refetch on curation events that recompute indices.
 */
export interface SseEvent {
  id: number;
  kind: string;
  entity_type: string;
  entity_id: number;
  actor?: string | null;
  client_id?: string | null;
  client_op_id?: string | null;
  ts?: string;
  payload?: unknown;
  post_state?: { analysis_inputs_hash: string; indices: unknown[] };
}

/**
 * A mutator describes one queue-able operation:
 * - `kind` discriminates against OpKind
 * - `onMutate` writes the optimistic cache effect; returns the rollback ctx
 * - `request` issues the HTTP call; honours AbortSignal
 * - `onSuccess` applies the server response to the cache
 * - `affectsExposurePeaks` (optional) tells hooks whether this op should
 *   register as a "peak op" against an exposure for `useExposureHasPendingPeakOps`
 */
export interface Mutator<TPayload, TResponse> {
  kind: OpKind;
  onMutate: (payload: TPayload, qc: QueryClient) => RollbackContext;
  request: (payload: TPayload, signal: AbortSignal) => Promise<TResponse>;
  onSuccess: (payload: TPayload, response: TResponse, qc: QueryClient) => void;
  affectsExposurePeaks?: (payload: TPayload, exposureId: number) => boolean;
}
