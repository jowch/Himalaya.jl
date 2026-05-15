import type { QueryClient } from "@tanstack/react-query";
import type { Comparison } from "../../api";

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
 * The set of mutation kinds the queue framework can carry.
 *
 * Most kinds are 1:1 with the server-side `user_actions.action` event kind
 * (e.g. `peak_added` → `apply_event!(..., kind="peak_added")`), but a few
 * are queue-side aliases for routes that emit a differently-named event:
 *
 * - `delete_index` → emits `speculative_deleted` (the only deletable index
 *   kind today is speculative; the OpKind name describes the *user gesture*,
 *   the event name describes the *durable mutation*).
 * - `reanalyze_exposure` → emits `analyze_run` (manual-triggered analyze).
 *
 * Kind-set membership predicates (`useExposureHasPendingPeakOps`,
 * `PEAK_AFFECTING_KINDS`) operate on OpKind, not the event name. If a new
 * mutator's user-gesture name diverges from its event kind, document the
 * mapping here and double-check the membership predicates.
 */
export type OpKind =
  | "peak_added" | "peak_excluded" | "peak_unexcluded" | "peak_removed"
  | "index_confirmed" | "index_unconfirmed"
  | "speculative_created" | "speculative_deleted"
  | "set_exposure_status" | "select_exposure"
  | "add_tag" | "remove_tag"
  | "post_message" | "update_sample"
  | "reanalyze_exposure"
  | "delete_index"
  // Compare page (Plan §Phase 3). One queue mutator handles both create
  // (`POST /api/comparisons`) and submit (`POST /api/comparisons/:id/submit`)
  // — `payload.id` discriminates. The wire-side event names
  // (`comparison_created`, `comparison_submitted`) diverge from the OpKind
  // because a single user gesture (Save) can map to either event.
  // `comparison_delete` maps 1:1 with `comparison_deleted`.
  | "comparison_save" | "comparison_delete";

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
  /** AbortController for the in-flight HTTP request; abort when SSE resolves first. */
  controller?: AbortController;
}

/**
 * Curation-frame `post_state`: the recomputed indices snapshot threaded onto
 * `peak_*` / `analyze_run` SSE frames so the cache can replay without a
 * refetch.
 */
export interface CurationPostState {
  analysis_inputs_hash: string;
  indices: unknown[];
}

/**
 * The shape of an SSE frame as parsed from the JSON `data:` line. Mirrors
 * the Julia-side `broadcast_event!` JSON shape (events.jl).
 *
 * Optional `post_state` carries an enriched snapshot for replay-without-
 * refetch. Its shape depends on `kind`: curation events (peak_*, analyze_run)
 * carry a `CurationPostState`; comparison events (comparison_created /
 * comparison_submitted) carry the full `Comparison` projection — the same
 * shape `fetch_comparison_with_members` / `GET /api/comparisons/:id` returns
 * (per Compare UX A-5 Step 5b). Consumers narrow by `kind` and cast.
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
  post_state?: CurationPostState | Comparison;
}

/**
 * The framework-owned fields that the replay coordinator threads into every
 * synthesized response. Mutators MAY destructure these off via
 * `stripQueueMetadata`. The values come from the SSE frame, not the payload.
 */
export interface QueueResponseMeta {
  event_id: number;
  client_op_id: string | null | undefined;
  analysis_inputs_hash: string | undefined;
  view_row_id?: number;
}

/**
 * The flat payload shape `useQueueMutation` actually constructs and hands to
 * mutator callbacks: the OpPayload metadata fields, the scope object the hook
 * caller passed in, and the per-call input merged together. Defined once so
 * mutators don't have to re-derive it via `as unknown as` casts.
 */
export type FlatPayload<TInput, TScope> = OpPayload<TInput> & TScope & TInput;

/**
 * A mutator describes one queue-able operation:
 * - `kind` discriminates against OpKind
 * - `onMutate` writes the optimistic cache effect; returns the rollback ctx
 * - `request` issues the HTTP call; honours AbortSignal
 * - `onSuccess` applies the server response to the cache
 * - `synthesizeFromSse` (optional) builds the SSE-wins synthetic response
 * - `affectsExposurePeaks` (optional) tells hooks whether this op should
 *   register as a "peak op" against an exposure for `useExposureHasPendingPeakOps`
 *
 * Three generics: `TInput` (per-call input), `TScope` (closure-injected at hook
 * construction), `TResponse` (server response). Callbacks receive the flat
 * merged payload — no casts at the consumer layer.
 */
export interface Mutator<TInput, TScope, TResponse> {
  kind: OpKind;
  onMutate: (payload: FlatPayload<TInput, TScope>, qc: QueryClient) => RollbackContext;
  request: (payload: FlatPayload<TInput, TScope>, signal: AbortSignal) => Promise<TResponse>;
  onSuccess: (payload: FlatPayload<TInput, TScope>, response: TResponse, qc: QueryClient) => void;
  /**
   * Build a synthetic response for the SSE-wins path. When the SSE for our
   * own pending op lands before the HTTP response, the deferred resolves
   * with this value so `onSuccess` can apply the cache effect without
   * waiting for the (now-aborted) HTTP call. Return `undefined` to fall
   * back to the generic `{ ...base, ...payload }` shape.
   *
   * The shape contract for what this returns IS the same shape `onSuccess`
   * expects in its second argument — co-locating the two prevents drift.
   */
  synthesizeFromSse?: (remote: SseEvent, base: QueueResponseMeta) => TResponse | undefined;
  affectsExposurePeaks?: (payload: FlatPayload<TInput, TScope>, exposureId: number) => boolean;
  /**
   * When set, an HTTP 404 response is treated as a no-op success: the
   * framework skips both the toast and the rollback, leaving the optimistic
   * cache state in place. Use on idempotent remove/delete operations where
   * "the server doesn't have it" is the desired end state.
   *
   * Motivation: under 5xx-then-retry, a successful first attempt deletes the
   * row; the client retries on the missed response and gets a 404 from a
   * server that's already done the work. Without this flag, the rollback
   * re-inserts the row that the server (and SSE) have already removed,
   * producing a phantom row visible until the next refetch.
   */
  treats404AsSuccess?: boolean;
}
