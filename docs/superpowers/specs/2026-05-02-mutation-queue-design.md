# Server-Reconciliation Mutation Queue — Design Spec

**Date:** 2026-05-02
**Status:** Draft

## Context

HimalayaUI today is correctness-first multiplayer with last-writer-wins resolution. Plan 7 (R0–R5a) shipped the foundation — auto/curation peak split, content-hash memoization, structured event log, SSE broadcast — but the user-perceived snappiness ceiling is set by two structural facts:

1. **The autoReanalyze double round-trip.** Every peak edit chains `peak op → reanalyze → invalidate` ([queries.ts:103-152](../../packages/HimalayaUI/frontend/src/queries.ts)). The user's click pays two HTTP round-trips before indices reflect their edit. The recently-added 2000ms debounce on `StaleIndicesBanner` ([commit 824a105](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx)) is a workaround for the resulting flicker — it masks symptoms, not the cause.

2. **TanStack Query's invalidate-and-refetch model assumes server-truth-as-displayed-state.** The cache lags every mutation by one round-trip. For a SAXS curation tool used in long sessions, this lag accumulates into perceived heaviness even on healthy networks.

Plan 7 R5b (`If-Match` headers + 409 retry on conflict) was deferred behind an instrumentation gate (≥2% delta-event collision rate over ≥4 weeks / ≥500 events). The collision-resolution problem R5b was meant to solve is real but small. The bigger UX win — eliminating the autoReanalyze double round-trip and making every interaction feel immediate — needs a different architecture, not a refinement of R5b.

This spec proposes that architecture: a **server-reconciliation mutation queue**. The pattern is well-trodden — Linear, Replicache, Convex, TanStack DB all instantiate it. Server holds authoritative state `S`. Client maintains a queue of pending semantic operations `L₁…Lₙ`. Rendered state = `S` + queue, applied optimistically and instantly on click. On incoming SSE event: undo queue → apply remote → re-run pending mutators. On HTTP confirmation: pop matching op from queue. On failure: pop and roll back.

The pattern fits because Himalaya's user operations are coarse semantic actions (`peak_added`, `index_confirmed`), not character-level edits. CRDT and OT machinery are unnecessary; the existing `apply_event!` dispatcher contract on the server already exhibits the same purity properties the client-side queue needs.

### Relationship to prior spec / R5b disposition

R5b's instrumentation gate was correct for *its* design — building 409+retry without evidence of conflict frequency would have been speculative. The queue model proposed here is a different architecture whose primary value is **elimination of the autoReanalyze double round-trip and structural correctness of optimistic UI**. Those benefits are independent of conflict frequency.

The queue model still resolves conflicts — via replay-as-rerun (mutators re-run against new base on remote events). For *additive* curation (two users adding different peaks, excluding different peaks), conflicts resolve invisibly: each mutator re-runs against the merged base and either applies or no-ops. For *toggle races* (A confirms an index while B unconfirms it), replay-rerun produces correct *final UI state* but the event log shows a flap that didn't reflect either user's discrete intent — observable in instrumentation, not in UX.

R5b is therefore not deferred but **redirected**: the conflict-resolution problem doesn't go away, it gets a different shape. Replay-rerun handles the additive case better than 409+retry would (no user-visible interruption); the toggle case is no worse and arguably more honest (the flap is auditable). The instrumentation that R4 ships (`analyze_run` payloads, conflict-relevant metrics) remains valuable for observability and surfaces toggle-flap rates if they matter.

This is a judgment call, not data-driven, and the spec should own that. See [Fallback triggers](#fallback-triggers) for the pre-commitments that make this judgment falsifiable in implementation.

## Goals

- Every user mutation appears applied to the UI within one frame of click. No "your edit failed, retry" friction outside genuine failure cases.
- The autoReanalyze double round-trip is structurally eliminated. Server reanalyzes idempotently inside the same handler that applies the curation; client makes one round-trip, gets post-analysis state in the response.
- Conflicts between concurrent users resolve invisibly when the resolution is unambiguous (most cases), and surface comprehensible affordances when not (replay alters my optimistic effect, ambient state changed under a hovered element).
- Network failure is graceful: optimistic UI persists across in-flight retries; user sees "saving…" indication only past a perceptual threshold; queue survives reload within a session.
- All 16 mutations route through one uniform machinery. No queue-vs-standard branching at consumer sites.
- The architecture is symmetric with `apply_event!`: the same `(kind, entity_id, payload, base_state) → effects` shape on both sides of the wire.

## Non-goals

- **Cross-session undo (Cmd+Z).** The primitives exist (compensating ops via `peak_unexcluded`, `index_unconfirmed`); the UI is its own design and benefits from but doesn't require this spec. Scoped out for clean focus; addressed in a follow-up.
- **`Last-Event-ID` server-side SSE replay.** With HTTP-authoritative confirmation (see Q3 below), SSE drops only affect cross-tab/user notification, which already reconciles via TanStack Query refetch on reconnect. Last-Event-ID is unnecessary for queue correctness.
- **Offline editing.** Persistence is sessionStorage-scoped; tab close drops pending ops. A real offline story would require localStorage persistence, schema-versioned op replay, and conflict resolution against significant server divergence — out of scope.
- **CRDT or OT machinery.** Operations are coarse semantic actions; the queue's replay-as-rerun is sufficient.
- **Automatic ML-driven conflict resolution.** Conflicts surface honest UX affordances; the system does not infer user intent beyond mutator rerun.
- **Replacement of TanStack Query.** The queue layers on top of `MutationCache` rather than competing with it.

## Pain points in the current design

Three load-bearing examples of the structural problem:

1. **Every peak edit is two round-trips.** [queries.ts:116-152](../../packages/HimalayaUI/frontend/src/queries.ts) — `useAddPeak`, `useRemovePeak`, `useSetPeakExcluded` all chain `peak op → autoReanalyze → invalidate`. The first round-trip writes the curation. The second triggers reanalysis. Indices don't update until both complete. Even on a 50ms-RTT network this is 100ms+ of latency on every peak edit, with two cache-invalidation cycles.

2. **The `StaleIndicesBanner` debounce is a workaround for flicker.** [StaleIndicesBanner.tsx:9](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx) introduces a 2000ms delay specifically to mask the transient mismatch between server-applied curation and not-yet-applied reanalysis. The banner is correct; the underlying flicker is the symptom of the deeper architectural mismatch.

3. **R5b would solve a narrower problem.** `If-Match` + 409 retry (~150 lines across ~8 routes) is a precise fix for the conflict-resolution case the prior spec gated behind instrumentation. It does nothing for the autoReanalyze double round-trip, the cache-flicker debounce workaround, or the broader "every interaction lags by one round-trip" feel. Investing the same engineering energy in the queue model yields all of these at once, with R5b's conflict-resolution problem subsumed as a side effect of replay-as-rerun.

A capability the current model can't reach:

4. **No queue-aware optimistic primitive exists.** TanStack Query's `onMutate` handles per-mutation optimistic updates in isolation, but there is no machinery for *queue-aware* optimistic UI: surviving remote events landing mid-flight, persisting across reload, applying uniform pending indication to all in-flight ops, and replaying mutators against new base state when the cache shifts under them. Building this on top of `onMutate` is the queue layer regardless of which architectural pattern we adopt.

## Design overview

Three architectural pieces, four implementation milestones.

### Architectural pieces

| Piece | Server-side | Client-side |
|---|---|---|
| **Identity** | `user_actions.client_op_id TEXT` | `crypto.randomUUID()` per op, `X-Client-Op-Id` header |
| **Idempotency** | Route-level wrapper + response-body cache (`idempotent_responses`); per-op-id lock serializes concurrent retries | Pop queue on first signal (HTTP or SSE); idempotent on second |
| **Confirmation** | HTTP response carries `{event_id, view_row_id, analysis_inputs_hash, peak_or_curation_row}` | Mutator's `onSuccess` updates cache from response |
| **Reconciliation** | SSE frame extended with `client_op_id`, `ts`, and `post_state` (on reanalysis events) | Replay coordinator: pop on `client_op_id` match, else rollback queue → apply remote → re-run mutators |
| **Persistence** | (existing) `user_actions` event log + `idempotent_responses` cache (TTL'd) | sessionStorage of `{schemaVersion, kind, payload, clientOpId}` |

### Architectural decisions (digest)

The fourteen open questions raised pre-spec resolved as follows:

| # | Question | Decision |
|---|---|---|
| Q1 | Unit of work | Semantic operations (not request descriptors) |
| Q2 | Replay semantics | Replay-as-rerun: mutator is `(op, baseState) → (cacheEffect, request)`, runs again on remote events |
| Q3 | Confirmation primitive | HTTP-authoritative; SSE secondary; `client_op_id` per request |
| Q4 | Library choice | Roll-our-own on TanStack Query `MutationCache` |
| Q5 | Server idempotency | Route-level wrapper; response-body cache (`idempotent_responses` table) keyed by `client_op_id`; per-op-id in-process locking serializes concurrent retries |
| Q6 | Persistence | sessionStorage; matches `client_id` lifetime |
| Q7 | Optimistic seam | Cache-merged via `setQueryData`; `MutationCache` is the queue |
| Q8 | Stale banner | Gate on `useExposureHasPendingPeakOps`; retain short (~150ms) debounce as backstop for replay-coordinator-induced cross-entity refetch races |
| Q9 | Failure taxonomy | Three classes: validation, conflict, infrastructure |
| Q10 | Queue scope | Uniform: all 16 mutations queue-shaped; null mutators for ops without optimistic effect |
| Q11 | Server-computed reads | Gate on `useExposureHasPendingPeakOps` |
| Q12 | Two-tab same-username | Resolved by PR #32; no further work |
| Q13 | R5b disposition | Redirected (not deferred): this spec subsumes R5b's conflict-resolution problem via replay-as-rerun, with different tradeoffs (see Context §Relationship to prior spec) |
| Q14 | SSE frame extension | `client_op_id` and `ts` added; `post_state` (containing post-reanalyze `analysis_inputs_hash` and updated indices array) added on curation events that triggered reanalysis |

### Implementation milestones

| Milestone | Independent value | Multiplayer enabler? |
|---|---|---|
| M0 | Schema + idempotency primitives + `analyze_exposure!` fast-skip refactor + response-body cache. | Mostly backend; minimal frontend (`clientOpId` mint helper + header threading). Behavioral no-op for users. |
| M1 | Queue infrastructure. Mutator framework with deferred-promise resolution, replay coordinator, sessionStorage persistence. | Frontend infrastructure; no behavior change yet. |
| M2 | Mutator migration. All 16 mutations migrated to the queue framework. autoReanalyze chain dissolved. Stale banner gating updated. Failure-class indicators (validation toast, infrastructure banner) shipped alongside. | Primary deliverable; user-visible weightless UX. |
| M3 | Cleanup. Drop unused fallback paths in PR #32's filter; consolidate hash memoization invariants. | Code health. |
| (Follow-up) | Polished UI affordances. Per-element pending state, aggregate "saving…" indicator, ambient-change highlight. | Separate PR; the *failure-class* indicators ship in M2, the *success-class* polish is the follow-up. |

Each milestone is mergeable independently. M0 ships behavioral no-op. M1 ships infrastructure with no consumer. M2 is where user-visible behavior changes. M3 is mechanical.

---

## M0: Schema, idempotency, and analyze fast-skip

M0 is larger than initially scoped. Three pieces all land before any frontend queue work:

1. **Schema and idempotency primitives** — `client_op_id` column, response-body cache, route-level wrapper.
2. **`analyze_exposure!` fast-skip refactor** — load-bearing for M2's snappiness claim. Without this, "synchronous reanalyze in the curation handler" relocates latency rather than eliminating it.
3. **`analyze_run` broadcast suppression** — when both skip flags fire, suppress the SSE broadcast (no observable state change → nothing to notify; avoids O(N²) frame fan-out under N concurrent tabs).

### Schema deltas

```sql
-- Add request-level idempotency identity to the event log.
-- Nullable for system events (pipeline, _system_request()) and pre-feature rows.
ALTER TABLE user_actions ADD COLUMN client_op_id TEXT;

-- Partial index: NULL is exempt; non-NULL gets indexed for lookup-by-op-id.
-- NOT UNIQUE: a single request may produce multiple events under one client_op_id
-- (speculative POST in routes_analysis.jl emits N events for N (phase, anchor)
-- combinations). Idempotency lives at the route handler level, not the row level.
CREATE INDEX IF NOT EXISTS idx_events_by_client_op_id
  ON user_actions(client_op_id) WHERE client_op_id IS NOT NULL;

-- Response-body cache for idempotent retries.
-- Why a separate table rather than reconstruction-from-event-log:
-- - Routes that emit N events under one request (speculative POST) would need
--   per-route reconstruction logic that reassembles array-shaped responses in
--   the original event order.
-- - Response shapes evolve; deployed clients retrying old client_op_ids
--   against an upgraded server need byte-identical replay of the original
--   response, not a fresh reconstruction.
-- - Storage cost is bounded: ~1KB per cached response × ops-per-session ×
--   active users. At lab scale this is ~MB, not a concern.
CREATE TABLE IF NOT EXISTS idempotent_responses (
    client_op_id  TEXT PRIMARY KEY,
    status_code   INTEGER NOT NULL,
    body          TEXT NOT NULL,           -- JSON
    created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- TTL-based GC: rows older than 1 hour are eligible for deletion.
-- 1 hour is generous for retry windows (typical retry budget is seconds, not
-- minutes). Run on each open_db, plus a periodic sweep — rare enough that a
-- simple "DELETE WHERE created_at < datetime('now','-1 hour')" suffices.
CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created
  ON idempotent_responses(created_at);
```

Add to `migrate_schema!`'s idempotent stmts array ([db.jl](../../packages/HimalayaUI/src/db.jl)) following the same shape as PR #32's `client_id` migration.

### `analyze_exposure!` fast-skip refactor

The current implementation in [pipeline.jl:597-608](../../packages/HimalayaUI/src/pipeline.jl) computes `hash_trace_file(dat_path)` (full file read + SHA-256) and `load_dat(dat_path)` *unconditionally before* the skip predicates fire. The skip flags only avoid CPU work for `findpeaks` and `indexpeaks`; file I/O happens every call. For a typical SAXS .dat file this is 1-5ms on local SSD plus a synchronous SSE broadcast — invisible to backend P95 metrics, fatal to perceived snappiness on the curation hot path.

**The refactor:** restructure the function so that DB-only signals can prove no-change before any disk I/O.

```julia
function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String;
                            trace_known_unchanged::Bool=false)
    # Caller hint: peak-curation routes know the trace file wasn't touched.
    # Pass trace_known_unchanged=true to skip hash_trace_file entirely. The
    # default (false) preserves the current behavior for callers who don't
    # know — reingest, scheduled scans, manual reanalyze.

    # DB-only reads — no file I/O yet.
    dat_path           = resolve_trace_path(db, exposure_id, analysis_dir)
    stored_trace_hash  = read_trace_hash(db, exposure_id)
    stored_inputs_hash = read_inputs_hash(db, exposure_id)
    autopeaks_count    = count_auto_peaks(db, exposure_id)
    indices_count      = count_indices(db, exposure_id)

    if trace_known_unchanged
        new_trace_hash = stored_trace_hash
        findpeaks_skipped = autopeaks_count > 0
    else
        new_trace_hash = hash_trace_file(dat_path)  # disk I/O
        findpeaks_skipped = (stored_trace_hash == new_trace_hash) && (autopeaks_count > 0)
    end

    # Hot-path optimization: if findpeaks is skipped AND we can compute the new
    # inputs_hash from auto_peaks + curations directly (no load_dat needed for
    # peak coordinates — `add` curations need sharpness from the trace, but only
    # when present), skip load_dat entirely.
    if findpeaks_skipped && !any_add_curations(db, exposure_id)
        new_inputs_hash = hash_peak_set_from_db(db, exposure_id)
        indexpeaks_skipped = (stored_inputs_hash == new_inputs_hash) && (indices_count > 0)
        if indexpeaks_skipped
            # Both stages skipped, no disk I/O. Fast path.
            # apply_event! is called normally; the broadcast-suppression
            # logic in events.jl detects the both-skipped case from the
            # payload and elides the SSE frame (see next subsection).
            apply_event!(db, _system_request();
                kind        = "analyze_run",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = Dict(
                    :findpeaks_skipped  => true,
                    :indexpeaks_skipped => true,
                    :trace_hash         => new_trace_hash,
                    :inputs_hash        => new_inputs_hash,
                    # ... post_state_size_bytes etc.
                ))
            return
        end
    end

    # Slow path: load trace, run findpeaks/indexpeaks as needed.
    q, I, σ = load_dat(dat_path)
    # ... existing code path emits its own analyze_run event with the
    # appropriate skip flags; broadcast fires normally for slow-path events.
end
```

Helpers:

- `read_trace_hash`, `read_inputs_hash`, `count_auto_peaks`, `count_indices`: small SELECTs returning scalars; trivial.
- `any_add_curations(db, exposure_id)`: `SELECT 1 FROM peak_curations WHERE exposure_id = ? AND kind = 'add' LIMIT 1`. Trivial.
- `hash_peak_set_from_db(db, exposure_id)`: reads `auto_peaks` (q, sharpness) joined against `peak_curations(kind='exclude')` for the exposure, computes `hash_peak_set` over the resulting (q, sharpness) tuples *without loading the trace*. Only valid when no `add` curations exist (those need sharpness from `Himalaya.sharpness(I)` — see [pipeline.jl:474](../../packages/HimalayaUI/src/pipeline.jl) and the `effective_peaks` contract in [event-log.md](../event-log.md) §1).

The `add`-curations case still needs `load_dat`. That's acceptable — the typical curation flow is exclude (>>far more common than add), and the spec explicitly accepts the slow path for the rarer add case. Future optimization: persist `add`-curation sharpness when `auto_peaks` are last computed and validate against trace_hash, eliminating the load_dat dependency entirely. Out of scope for M0.

### `analyze_run` broadcast suppression

```julia
# events.jl
if isdefined(@__MODULE__, :broadcast_event!)
    # Suppress broadcast for analyze_run events that recorded no observable
    # state change (both findpeaks and indexpeaks skipped). Other clients have
    # nothing to refetch; broadcasting just adds noise and triggers spurious
    # cache work in their TanStack Query layers.
    suppress = kind == "analyze_run" &&
               payload !== nothing &&
               get(payload, :findpeaks_skipped, false) &&
               get(payload, :indexpeaks_skipped, false)
    if !suppress
        broadcast_event!(...)
    end
end
```

Without this, a peak edit fires both a curation event (peak_added/peak_excluded/...) AND an analyze_run event, both broadcast to all subscribers. Under N concurrent tabs editing, broadcast volume is O(N²). The suppression reduces steady-state broadcast volume to one frame per user-meaningful event.

### `apply_event!` extension

Mirror PR #32's `client_id` plumbing exactly:

```julia
# events.jl — add client_op_id parameter alongside the existing client_id
function apply_event!(db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing)
    username     = get_username(req)
    client_id    = get_client_id(req)
    client_op_id = get_client_op_id(req)  # NEW: extracts X-Client-Op-Id header
    user_id      = username === nothing ? nothing : get_or_create_user!(db, username)
    # ... INSERT writes client_op_id alongside client_id
end
```

`get_client_op_id` is a 3-line helper in [actions.jl](../../packages/HimalayaUI/src/actions.jl) modeled on `get_client_id` from PR #32.

### Route-level idempotency wrapper

The wrapper records the route's *response body* keyed by `client_op_id` on first execution, and returns the cached body verbatim on retry. This is robust to multi-event routes (speculative POST emits N events under one `client_op_id`) and to response-shape evolution (a deployed-then-upgraded server returns byte-identical original responses to retry attempts using old `client_op_id`s).

```julia
# events.jl — new wrapper.

# Per-op-id locks serialize concurrent retries within the process. Two retries
# for the same client_op_id must not both run f() — the body emits events and
# we'd get duplicate downstream effects (two peak_curations rows for one logical
# operation). HimalayaUI is single-process per CLAUDE.md "one DB" model, so an
# in-memory ReentrantLock per op_id is sufficient; cross-process protection is
# not in scope.
const OP_LOCKS    = Dict{String, ReentrantLock}()
const OP_LOCKS_MU = ReentrantLock()

function _op_lock(op_id::String)::ReentrantLock
    lock(OP_LOCKS_MU) do
        get!(OP_LOCKS, op_id, ReentrantLock())
    end
end

"""
    with_idempotency(db, req) do
        # body runs at most once per client_op_id; otherwise the cached
        # response is returned verbatim.
    end

Request-level idempotency via response-body cache + per-op-id locking.

Behavior:
- NULL or missing `X-Client-Op-Id` → passthrough; body always executes,
  nothing cached. Preserves unit tests, `_system_request()` flows, and
  pre-feature client compatibility.
- Cache hit → return the previously-recorded `(status_code, body)` exactly,
  with no DB writes.
- Cache miss → acquire per-op-id lock, re-check cache (another request may
  have raced ahead and written the response), then execute the body and
  on success (status_code < 400) write to `idempotent_responses`.

The successful-only caching policy means: if the body's `SQLite.transaction`
rolls back (constraint violation, exception), no row in `user_actions` is
written and no row in `idempotent_responses` is written. Retry naturally
re-executes against fresh state. A failed first attempt + retry produces one
successful event, not zero or two.

If the transaction succeeds but the route handler errors *after*
apply_event! returns (formatting bug, downstream issue), the event is durable
in `user_actions` but the response cache row is not written. Retry will
re-execute the body — which on the second attempt sees the prior event row
and (depending on route logic) may either replay successfully or detect
a duplicate. This rare edge case warrants a route-level test per route.
"""
function with_idempotency(f, db::SQLite.DB, req::HTTP.Request)
    op_id = get_client_op_id(req)
    op_id === nothing && return f()

    # First fast-path check (lock-free read).
    cached = _lookup_cached_response(db, op_id)
    cached === nothing || return cached

    # Acquire per-op-id lock; re-check inside the lock (another request may
    # have written between our check and our acquisition).
    lock(_op_lock(op_id)) do
        cached2 = _lookup_cached_response(db, op_id)
        cached2 === nothing || return cached2

        response = f()
        if response.status < 400
            DBInterface.execute(db,
                "INSERT INTO idempotent_responses (client_op_id, status_code, body) VALUES (?, ?, ?)",
                [op_id, Int(response.status), String(response.body)])
        end
        return response
    end
end
```

**Cache-coherence note for cached responses.** The cached body is the response that the *original* request produced. If, between the first commit and a subsequent retry, another user deletes an entity referenced in that body, the cached response still mentions the now-deleted entity. The retry's `onSuccess` writes that body into the client's cache, then the SSE event for the deletion arrives shortly after (already enqueued at the time of deletion) and reconciles the cache. There is a brief window where the client's cache is stale — bounded by SSE delivery (typically <1s) — and self-heals via the replay coordinator. We document this rather than try to invalidate the response cache on entity deletion: building a per-entity reverse index from `client_op_id` to referenced entities would be substantial machinery for a self-healing edge case.

For the **rehydrate** path specifically (sessionStorage queue replayed on tab reawake): the freshly-fetched cache on rehydrate is the source of truth; the cached response from a retry is layered on top. If the cached response references an entity that's been deleted, the brief window of staleness is the same as the in-session retry case and self-heals identically. Rehydrate intentionally relies on HTTP-only confirmation (the SSE frames from the original session were missed during sleep, by design — see Q6).

**Why response caching as the primary path, not reconstruction.** An earlier draft of this spec proposed reconstructing the response from the event log on retry. This works for routes that emit a single event with a simple response shape, but fails on:

- **Multi-event routes.** [routes_analysis.jl:264-310](../../packages/HimalayaUI/src/routes_analysis.jl) emits N events for N (phase, anchor) combinations chosen by the client. Reconstruction would need to assemble an array-shaped response in original event order, traversing `user_actions` joined to `view_row_id` lookups in `index_groups`/`index_group_members`/`indices`.
- **Response-shape evolution.** A client retrying an old `client_op_id` against an upgraded server expects the response shape it would have received at original-call time, not the new shape. Reconstruction always synthesizes against current schema.
- **Per-route maintenance burden.** Reconstruction code is structurally similar to the cache-invalidation code we're trying to escape.

Response-body caching costs ~1KB per cached response × ops-per-session × active users. At lab scale (single-digit users, hundreds of ops per session) this is single-digit MB. TTL'd to 1 hour, never grows unbounded. The per-op-id `OP_LOCKS` map grows similarly bounded; lazy GC alongside the response cache TTL.

### Frontend: client_op_id minting

```typescript
// lib/clientOpId.ts — new file
export function newClientOpId(): string {
  return crypto.randomUUID();
}

// api.ts — extend AuthOpts and request()
export interface AuthOpts {
  username?: string;
  clientId?: string;
  clientOpId?: string;  // NEW
}

// In request(), threads X-Client-Op-Id on mutations only (matches X-Username, X-Client-Id rule)
```

### Wire format

| Header | Source | Lifetime | Purpose |
|---|---|---|---|
| `X-Username` | persisted Zustand state | session login → logout | audit attribution |
| `X-Client-Id` | sessionStorage UUID (PR #32) | tab session | SSE routing/filtering |
| `X-Client-Op-Id` | per-mutation UUID | one request | request-level idempotency, queue confirmation |

Headers compose. A single mutation request carries all three.

### SSE frame extension

```json
{
  "id":           1234,
  "kind":         "peak_added",
  "entity_type":  "exposure",
  "entity_id":    42,
  "actor":        "alice",
  "client_id":    "uuid-tab-A",
  "client_op_id": "uuid-op-789",
  "payload":      { "q": 1.234 },
  "ts":           "2026-05-02T19:23:45.123Z",
  "post_state": {
    "analysis_inputs_hash": "abc123...",
    "indices":              [ /* updated indices array */ ]
  }
}
```

`client_op_id` lets the queue's replay coordinator pop matching ops on SSE arrival. `ts` enables future audit/UI affordances. Both NULLable for pre-feature events and `_system_request()` flows.

`post_state` is the load-bearing addition for replay-without-refetch (see M1 replay coordinator below). Curation events that triggered a synchronous reanalysis in their handler include the resulting `analysis_inputs_hash` and the updated indices list. The replay coordinator applies these to the cache atomically — no separate refetch round-trip, no cross-entity refetch race, no banner flicker. Events that didn't trigger reanalysis (tag adds, chat messages) omit `post_state`. The field is omitted, not null, so existing serialization code stays minimal.

Frame size grows by ~1-3KB per curation event for typical exposures, but the indices array is unbounded — a heavily-curated exposure with many speculative-builder anchors can drive `post_state` past 8KB. Under N concurrent tabs the broadcast volume is N × frame size, so the gain from M0's `analyze_run` broadcast suppression can be partly clawed back by `post_state` enlargement.

To make this measurable rather than speculative: **M0 includes payload-size observability.** `analyze_run` event payloads gain a `post_state_size_bytes` field (computed cheaply from the JSON-serialized indices array length). M2 reads these to decide whether `post_state` is sustainable as designed or needs the fallback in Open Question #4 (compact "you should refetch the following keys" payload, accepting the brief banner-gating window).

### M0 test impact

- Julia: tests for `client_op_id` migration, `get_client_op_id` extraction, `apply_event!` persistence (header present + absent), `with_idempotency` cache-hit and cache-miss paths, failed-first-attempt-then-retry produces one successful event (not zero, not two), multi-event route (speculative POST) cache-hit returns identical body.
- Julia: tests for `analyze_exposure!` fast-skip — verify zero file I/O when `trace_known_unchanged=true` and no add-curations, verify slow path still correct when add-curations exist or trace changed.
- Julia: test for `analyze_run` broadcast suppression — both-skipped events don't reach SSE subscribers, single-skipped events still broadcast.
- Frontend: `clientOpId` mint helper, header threading in `request()`.
- Existing tests stay green; M0 is additive.

### M0 independent value

Mostly backend (the minimal frontend work is `clientOpId` minting + header threading; users see no behavioral change). After M0 ships:
- Every route is idempotent under retry.
- `analyze_exposure!` is genuinely cheap on the no-change path (microseconds, not milliseconds; no file I/O).
- `analyze_run` broadcasts only when something changed, eliminating the noisiest source of cross-tab spam.

These are independently valuable even if M1/M2 never ship — for example, the autoReanalyze chain in [queries.ts:103-114](../../packages/HimalayaUI/frontend/src/queries.ts) immediately benefits from the fast-skip refactor (each chain's second round-trip becomes effectively free when the curation didn't actually change the effective peak set, which is the idempotent re-exclude case).

---

## M1: Queue infrastructure

### Mutator interface

```typescript
// lib/queue/types.ts
export interface OpKind {
  // Op kinds correspond 1:1 with user_actions.action values for queue-shaped ops.
  // 'reanalyze_exposure' and 'select_exposure' have null mutators (Q10).
}

export interface Mutator<TPayload, TResponse> {
  kind: OpKind;
  // Optimistic side: applies the predicted effect to QueryClient cache.
  // Returns rollback context (the prior cache snapshot) for replay-rollback.
  onMutate: (
    payload: TPayload,
    qc: QueryClient,
  ) => RollbackContext;
  // Request side: builds the HTTP call.
  request: (payload: TPayload) => Promise<TResponse>;
  // Success side: writes server-truth response to cache, replacing the optimistic.
  onSuccess: (
    payload: TPayload,
    response: TResponse,
    qc: QueryClient,
  ) => void;
  // Predicate: does this op affect <exposureId>'s peak set?
  // Used by useExposureHasPendingPeakOps for Q11 gating.
  affectsExposurePeaks?: (payload: TPayload, exposureId: number) => boolean;
}
```

The mutator framework is a thin wrapper over TanStack Query's `useMutation`. The custom layer is:
- `client_op_id` minting at queue time
- Mutation-cache integration so replay coordinator can find pending ops
- sessionStorage persistence
- Replay coordinator (separate concern below)

### Replay coordinator

The coordinator uses `mutationCache.getAll()` directly — no parallel queue data structure. **One non-obvious dependency:** TanStack Query's `MutationCache` stores mutations in a `Set<Mutation>` and `getAll()` returns `Array.from(this.#mutations)`. JS `Set` preserves insertion order per ECMAScript spec, so iteration is in queue order. This isn't loudly documented but is a stable property of the implementation. The spec depends on it; the M1 test suite asserts it explicitly so we'd notice if a future TanStack version changed.

**Mutator confirmation cannot use `Mutation.continue()`.** TanStack Query v5's `.continue()` is for resuming paused mutations and re-runs `mutationFn` — it cannot mark a mutation as having succeeded out-of-band. There is no public API for "resolve this mutation as if it returned X." The canonical workaround is the **deferred-promise pattern**: `mutationFn` awaits a deferred promise that resolves on whichever signal arrives first (HTTP response or matching SSE event).

```typescript
// lib/queue/types.ts
interface PendingDeferred<T> {
  promise: Promise<T>;
  resolve: (value: T) => void;
  reject: (reason: unknown) => void;
}

// Module-scoped registry, keyed by client_op_id.
const pendingDeferreds = new Map<string, PendingDeferred<unknown>>();

function makeDeferred<T>(clientOpId: string): PendingDeferred<T> {
  let resolve!: (v: T) => void;
  let reject!: (e: unknown) => void;
  const promise = new Promise<T>((res, rej) => { resolve = res; reject = rej; });
  const d = { promise, resolve, reject };
  pendingDeferreds.set(clientOpId, d as PendingDeferred<unknown>);
  return d;
}

// In each mutator's mutationFn:
async function peakAddMutationFn(
  payload: PeakAddPayload,
  ctx: { signal: AbortSignal },
): Promise<PeakAddResponse> {
  const clientOpId = payload.clientOpId;
  const deferred = makeDeferred<PeakAddResponse>(clientOpId);
  // Race: HTTP response wins → resolve with body. Or replay coordinator
  // resolves via SSE match → resolve with synthesized response. Or HTTP
  // errors → reject.
  api.addPeak(...).then(response => deferred.resolve(response))
                  .catch(err => deferred.reject(err));
  // AbortSignal handling: if TanStack Query cancels the mutation (e.g.,
  // component unmount, tab close hook), reject the deferred so the registry
  // entry doesn't leak. Without this, an aborted fetch's .then/.catch may
  // never fire and the Map entry persists for the page lifetime.
  ctx.signal.addEventListener("abort", () => {
    deferred.reject(new DOMException("aborted", "AbortError"));
  });
  try {
    return await deferred.promise;
  } finally {
    pendingDeferreds.delete(clientOpId);
  }
}
```

```typescript
// lib/queue/replayCoordinator.ts
function handleRemoteEvent(remote: SseEvent, qc: QueryClient, mc: MutationCache) {
  // SSE-confirms-my-pending-op path: resolve the deferred. mutationFn awakens,
  // returns to TanStack Query, which fires onSuccess and pops the mutation.
  if (remote.client_op_id) {
    const deferred = pendingDeferreds.get(remote.client_op_id);
    if (deferred) {
      deferred.resolve(synthesizeResponseFromSse(remote));
      return;
    }
  }

  // Remote event is from another tab/user. Walk pending mutations.
  const pending = mc.getAll().filter(m => m.state.status === "pending");

  // Roll back optimistic effects in reverse queue order. Each mutator's
  // onMutate stored a rollback context (the prior cache snapshot); applying
  // it restores baseline.
  for (const m of [...pending].reverse()) {
    const ctx = m.state.context as RollbackContext | undefined;
    ctx?.restore?.();
  }

  // Apply remote event to cache. For curation events the SSE payload
  // carries enough to update peaks/curations; for analyze_run events the
  // payload carries the new analysis_inputs_hash so we can update both
  // exposure-entity and indices-list cache atomically (no separate refetch
  // round-trip, avoiding the cross-entity refetch race that would otherwise
  // surface a banner flicker — see Q8).
  applyRemoteToCache(remote, qc);

  // Re-run optimistic effects in queue order against new base.
  for (const m of pending) {
    m.options.onMutate?.(m.state.variables);
  }
}
```

**Why HTTP-and-SSE both feed the same deferred:** the spec settled in Q3 that either signal can confirm an op, idempotently. Concretely the deferred is resolved on first arrival; the second is a no-op (resolving an already-resolved promise has no effect).

**SSE payload enrichment for replay-without-refetch.** The current SSE frame for curation events carries the event itself (kind, payload). For the replay coordinator to update the cache without invalidating-and-refetching, the frame needs enough payload to apply the cache mutation directly. For peak ops this is the curation row + the new `analysis_inputs_hash` from the synchronous reanalyze. We extend the broadcast to carry these — the same shape that an `analyze_run` event would carry, but inlined into the curation event so clients have one atomic frame to apply rather than two ordered frames. This eliminates a class of cross-entity refetch race that would otherwise re-introduce banner flicker (see Q8).

### Persistence layer

```typescript
// lib/queue/persistence.ts
const STORAGE_KEY = "himalaya-ui:queue";
const SCHEMA_VERSION = 1;

interface PersistedOp {
  schemaVersion: number;
  kind: OpKind;
  payload: unknown;
  clientOpId: string;
}

// Hook into MutationCache subscribe() to mirror pending ops to sessionStorage.
// On QueryClient init: rehydrate by re-running each persisted op's mutator,
// retried with the original clientOpId. Server-side idempotency (M0) returns
// the cached response if the op already landed; otherwise applies fresh.
```

Rehydrate flow:

1. Reload TanStack Query cache from server (normal startup).
2. For each persisted op, run its mutator's `onMutate` against the freshly-fetched cache.
3. Re-fire the request with the original `client_op_id`. Server returns either the original response (idempotent hit) or applies the op fresh.
4. On all confirmations, sessionStorage clears.

Ops with mismatched `schemaVersion` are dropped with a one-time toast: "your edits from before the last update were lost."

**Rehydrate intentionally relies on HTTP-only confirmation.** The dual-signal race (HTTP-or-SSE-whichever-arrives-first) is for in-session ops only. On rehydrate, any SSE frames from the original session were missed during sleep — by design, since the spec opted out of `Last-Event-ID` server-side replay. The replay coordinator's SSE-confirms-pending-op codepath therefore doesn't trigger for rehydrated ops; only HTTP can confirm them. This is fine because (a) every queue-shaped mutation is idempotent on `client_op_id`, and (b) post-rehydrate SSE frames flow normally and reconcile any divergence the rehydrate response introduced.

### `useExposureHasPendingPeakOps` hook

A small composable used by both Q8 (stale-banner gating) and Q11 (speculative-snap query gating). One source of truth for "does this exposure have pending writes that affect its peak set?"

```typescript
// lib/queue/hooks.ts
const PEAK_AFFECTING_KINDS = new Set<OpKind>([
  "peak_added", "peak_excluded", "peak_unexcluded", "peak_removed",
  "reanalyze_exposure",
]);

export function useExposureHasPendingPeakOps(exposureId: number | undefined): boolean {
  const pending = useMutationState({
    filters: {
      predicate: (m) => {
        if (m.state.status !== "pending") return false;
        const op = (m.state.variables ?? m.options.meta) as OpPayload | undefined;
        if (op === undefined) return false;
        if (!PEAK_AFFECTING_KINDS.has(op.kind)) return false;
        return op.exposureId === exposureId;
      },
    },
  });
  return pending.length > 0;
}
```

The hook re-runs when `mutationCache` state transitions, which is fine — TanStack Query already optimizes this internally and consumers see only the boolean. Same primitive backs the deferred indicator follow-up's per-element pending state.

### Failure-mode taxonomy

Three classes, three different UX affordances. Implementation handles each by inspecting `mutation.state` and `mutation.failureReason`:

| Class | Server signal | Queue behavior | UX |
|---|---|---|---|
| **Validation** | HTTP 4xx (400 invalid input, 404 entity missing, 422 semantic rejection) | Pop, rollback, surface error | Toast + inline error |
| **Conflict** | Replay-rerun produced no-op or different effect | Pop silently if no-op; soft annotation if changed | Often invisible; subtle inline notice when relevant |
| **Infrastructure** | HTTP timeout / 5xx / network error | Retry with backoff (idempotent on `client_op_id`) | "Saving…" indicator past 500ms threshold; "couldn't save" past 30s |

The "ambient changed under you" case (remote event landed, my queue is empty) isn't a failure mode — it's normal cache update.

**Entity-no-longer-exists during persistence rehydrate.** If sessionStorage holds a queued op against an entity that's been deleted before the tab reawakened, the rehydrated request hits HTTP 404. This is **Validation** (the op is invalid in the current world), with a class-specific message: "the peak you tried to edit was removed by another user; your edit didn't apply." No new failure class needed — the same toast pattern applies, just with rehydrate-aware copy.

**Cascading rejections in the queue.** If a queued op fails Validation because its target entity was deleted, *other* queue ops that referenced the same entity (e.g., `index_confirmed` for an index whose anchor peak was the deleted one) will fail their own Validation when they reach the server. The replay coordinator's mutator-rerun against the post-deletion base correctly produces no-ops or new validation failures for these dependent ops; each failure surfaces as its own Validation toast. There is no "first failure cascades and aborts the rest of the queue" semantic — every op confronts the new world independently and reports its own outcome. This keeps the failure model legible (each toast names exactly one op) at the cost of toast volume during pathological multi-op-against-deleted-entity cases. Acceptable in practice; the alternative ("one cascade summary") is harder to make accurate.

**Indicators for failure classes ship in M2, not the follow-up.** Without at least the Validation toast and Infrastructure banner, network failures during the optimistic window are invisible to the user — the cache shows the optimistic state, the request fails silently, the user thinks they saved. The "weightless" goal in §Goals is undefined for this case unless the user has a way to know something failed. The polished success-state indicators (per-element pending pulse, aggregate "saving…" widget, ambient-change highlight) remain in the follow-up. The failure-class indicators are scoped into M2 as a load-bearing UX commitment.

### M1 test impact

- Frontend: queue framework unit tests (mutator interface, replay coordinator with synthetic SSE events, persistence rehydrate cycle).
- M1 ships no consumer migrations, so existing tests stay green.

### M1 independent value

Frontend infrastructure with no consumer. Validates the framework against synthetic ops before any real mutation depends on it. Same shape as the original Plan 7 R4 dispatcher infrastructure landing before consumers in Plan 7 R5a.

---

## M2: Mutator migration

All 16 mutations from [queries.ts](../../packages/HimalayaUI/frontend/src/queries.ts) migrate to the queue framework. The migration is per-mutation, but the order is determined by mutator complexity:

### Migration order and atomicity

The slices have different atomicity requirements. Some are independent; the peak-op slice is not.

1. **Trivial mutators first** — cache-write only, no inter-op interaction: tag add/remove (4 mutations), chat message post, sample update, exposure status, `useSelectExposure`. Each is independently mergeable. These exercise the framework on simple cases before the load-bearing peak-op work.

2. **Peak-op slice (one atomic PR).** Three mutations migrate together with their server-side prerequisites: `useAddPeak`, `useRemovePeak`, `useSetPeakExcluded`. The PR includes:
   - Backend: extend curation route handlers in [routes_peaks.jl](../../packages/HimalayaUI/src/routes_peaks.jl) to call `analyze_exposure!(db, id; trace_known_unchanged=true)` synchronously inside `with_idempotency`. Response shape extends to `{event_id, view_row_id, analysis_inputs_hash, peak_or_curation_row}`.
   - Backend: SSE broadcast for these events extends to carry the post-reanalyze `analysis_inputs_hash` and updated indices payload (per the replay-without-refetch requirement above).
   - Frontend: migrate the three mutators to the queue framework with deferred-promise resolution.
   - Frontend: delete `autoReanalyze` from [queries.ts:103-114](../../packages/HimalayaUI/frontend/src/queries.ts).
   - Frontend: update [StaleIndicesBanner.tsx](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx) to gate on `useExposureHasPendingPeakOps`; reduce debounce to ~150ms (animation-frame-class backstop).

   The peak-op slice **cannot be subdivided** without producing intermediate states where some mutations bypass the queue while others use it, breaking replay-coordinator invariants and producing user-visible flicker. The work is naturally one PR.

3. **Index/group ops** — independent: `useAddIndexToGroup`, `useRemoveIndexFromGroup`, `useDeleteIndex`.

4. **Speculative ops** — independent: `useCreateSpeculative`, `useDeleteIndex`. The speculative-snap query gates on `useExposureHasPendingPeakOps` per Q11. Specifies UI behavior during gate: query stays disabled (no fetch); the SpeculativeBuilder modal shows the **last successful response** with a small "updating to latest…" subtext until pending peak ops settle, then re-enables the query. Avoids the gate producing a blank or skeleton state for the brief window.

5. **Null/trivial-mutator ops** — independent: `useReanalyzeExposure` (null optimistic effect; see ordering note below).

### `useReanalyzeExposure` ordering with pending peak ops

Open question resolved here so the migration order is clear. When the user clicks "Re-analyze" while peak ops are pending in the queue, the explicit reanalyze enters the queue **after** them in insertion order. It does not short-circuit pending peak ops. Reasoning:

- Short-circuiting would mean reanalyze runs *before* pending peak ops apply server-side, producing analysis state that ignores them — the opposite of what the user is asking for.
- FIFO ordering means the reanalyze waits ~one round-trip per pending peak op (typically zero pending ops; rarely more than one).
- Under R3 fast-skip (M0), the reanalyze itself is a near-no-op when peak-op-induced reanalyses already produced the latest hash. So the queue serializes correctly without paying multiple full-analysis costs.

The button's loading state derives from `useMutationState` filtered to this op kind: `isPending` becomes `pendingReanalyzeMutations.length > 0`.

### Mid-migration risk during M2

Between the trivial-mutator PRs and the peak-op slice, the cache is half-optimistic, half-invalidate-and-refetch. The replay coordinator only sees pending ops for migrated mutations. A non-migrated mutation in flight when a remote event arrives gets the legacy invalidate path.

The risk surfaces if a non-migrated mutation has *latency* concerns under remote-event-during-flight. The trivial mutators (tags, messages, sample fields) don't — they're scalar/additive and the legacy refetch handles them correctly. So in practice the half-migrated state is benign for the trivial slice.

The peak-op slice is the load-bearing case, and it ships atomically, so no half-migrated state exists for it.

### Server-side route handler shape (peak-op slice)

For reference, the route shape that lands in the peak-op slice (Migration step 2):

```julia
# Before:
@post "/api/exposures/{id}/peaks" function(req, id::Int)
    apply_event!(...)  # peak_added
    json_response(peak)
end

# After:
@post "/api/exposures/{id}/peaks" function(req, id::Int)
    with_idempotency(db, req) do
        SQLite.transaction(db) do
            apply_event!(...)  # peak_added
            analyze_exposure!(db, id, analysis_dir;
                              trace_known_unchanged=true)  # M0 fast-skip
        end
        json_response(Dict(
            :peak                  => peak_row,
            :event_id              => event_id,
            :view_row_id           => view_row_id,
            :analysis_inputs_hash  => current_inputs_hash(db, id),
        ))
    end
end
```

`trace_known_unchanged=true` is safe for any curation route — these handlers don't touch the .dat file. The fast-skip refactor (M0) makes the synchronous call cheap on no-input-change cases (target: <100µs total wall-clock; falsified by the corresponding fallback trigger). Routes that don't affect the effective peak set (tags, messages, status) omit the `analyze_exposure!` call entirely; their HTTP response carries no `analysis_inputs_hash`.

### Cleanup steps inside M2 (within respective slices)

- Delete `autoReanalyze` from [queries.ts:103-114](../../packages/HimalayaUI/frontend/src/queries.ts) — lands with the peak-op slice.
- Update [StaleIndicesBanner.tsx](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx): gate on `useExposureHasPendingPeakOps`; reduce debounce from 2000ms to ~150ms — lands with the peak-op slice.
- Update [sseSubscriber.ts](../../packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts): replace the PR #32 `client_id`-based skip filter with the queue's `client_op_id`-based confirmation match. The `client_id` filter becomes the *fallback* path for non-queue echoes (system events, pre-feature events) — lands with M1 if the framework needs it, otherwise the trivial-mutator slice.

### M2 test impact

- Vitest: every mutator round-trip (optimistic effect → request → confirmation → cache settled state). Existing testsets that asserted "after mutate + invalidate + refetch, the cache contains X" rewrite to "after mutate, the cache contains X immediately (optimistic), and after settle, the cache contains X (server-truth)." The Vitest harness already supports both via `await waitFor(...)`.
- Vitest: replay coordinator handles synthetic remote events that touch each pending op kind.
- Vitest: persistence rehydrate covers all op kinds.
- Playwright: two-context multiplayer test — User A and User B both add peaks to the same exposure; both succeed (deltas don't overlap). User A and User B both exclude the same peak; both succeed (idempotent under hashing). Same shape as the prior spec's R5 test, now exercising replay-rerun rather than 409+retry.
- Julia: server-side test that `analyze_exposure!` is called synchronously inside curation routes with `trace_known_unchanged=true`; that the M0 fast-skip path actually fires when expected.

### M2 independent value

User-visible weightless UX. Every interaction feels instant. The autoReanalyze double round-trip is gone. The stale banner only appears for genuine staleness.

---

## M3: Cleanup

Mechanical work that lands after M2 stabilizes:

- Remove the now-redundant `client_id`-as-skip-filter codepath in `sseSubscriber.ts` if M2 fully migrated all mutation echoes to queue confirmation. (May retain as defense-in-depth — decide during implementation.)
- Audit remaining `qc.invalidateQueries` calls in mutation `onSuccess` handlers. Most should be replaced with direct `setQueryData` from the response (the cache-merged optimistic pattern). Invalidation remains correct for SSE-triggered cache updates.
- Update [docs/event-log.md](../event-log.md) to document the queue's relationship to the dispatcher contract.
- Add a new section to CLAUDE.md under "HimalayaUI gotchas" naming the queue invariants and the `useExposureHasPendingPeakOps` pattern.

---

## Schema deltas (consolidated)

```sql
-- M0: idempotency identity on the event log
ALTER TABLE user_actions ADD COLUMN client_op_id TEXT;
CREATE INDEX IF NOT EXISTS idx_events_by_client_op_id
  ON user_actions(client_op_id) WHERE client_op_id IS NOT NULL;

-- M0: response-body cache for idempotent retries
CREATE TABLE IF NOT EXISTS idempotent_responses (
    client_op_id  TEXT PRIMARY KEY,
    status_code   INTEGER NOT NULL,
    body          TEXT NOT NULL,
    created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_idempotent_responses_created
  ON idempotent_responses(created_at);
```

The queue's *frontend* persistence is sessionStorage, not server-side; the backend persistence above exists solely to make HTTP-level retries idempotent.

## API surface changes

- Every mutation route gains optional `X-Client-Op-Id` header handling via `with_idempotency`.
- Curation-affecting routes ([routes_peaks.jl](../../packages/HimalayaUI/src/routes_peaks.jl), [routes_analysis.jl](../../packages/HimalayaUI/src/routes_analysis.jl)) call `analyze_exposure!` synchronously and return post-analysis state.
- Response shape extended to include `analysis_inputs_hash` on routes whose mutations affect the effective peak set.
- SSE frame extended with `client_op_id`, `ts`, and (on curation events that triggered reanalysis) `post_state`.

No breaking changes for clients that don't opt into the queue: missing `X-Client-Op-Id` makes `with_idempotency` a passthrough; the SSE frame's new fields are additive.

## Migration plan

Each milestone is its own PR, mergeable independently, behavior-stable in isolation:

1. **M0** — schema column + `idempotent_responses` table, idempotency wrapper, `analyze_exposure!` fast-skip refactor, `analyze_run` broadcast suppression, `client_op_id` plumbing through `apply_event!` and SSE frame. Behavioral no-op for unchanged clients.
2. **M1** — queue framework module, replay coordinator (with deferred-promise resolution), persistence layer. No consumers; module is dead code until M2 starts.
3. **M2** — phased mutator migrations. Trivial first (tags, chat, sample updates), then **the peak-op slice as one atomic PR** (see §M2 step 2 for the explicit "cannot be subdivided" rationale), then index/group/speculative, then null mutators. Failure-class indicators ship alongside the peak-op slice.
4. **M3** — mechanical cleanup.

The on-disk DB heals on next `open_db` for any subset merged. M0's column and table are additive; pre-feature events stay valid.

## Open questions

Implementation-level questions to resolve during the planning phase. The architectural decisions above don't block on these.

1. **Schema-versioning policy for sessionStorage ops.** The simple choice is "drop on schema_version mismatch with a one-time toast." A more sophisticated approach migrates op shapes across versions. For sessionStorage scope, the simple choice is probably right — the persistence window is one continuous session, and schema changes happen between sessions. Decide in M1.

2. **Aggregate "saving…" indicator threshold.** 500ms is the working assumption; literature suggests 200-1000ms range. Tune in the polish-indicator follow-up PR with real network observations. (The Validation toast and Infrastructure banner ship in M2; the *aggregate* indicator is the polish piece.)

3. **Replay rollback context shape.** TanStack Query's `onMutate` returns a context value used by `onError` to roll back. For replay-rollback we need this same context plus pending-op identity. Whether to extend the existing context shape or carry rollback state on `mutation.options.meta` is a small implementation detail. Decide in M1.

4. **SSE payload shape for replay-without-refetch.** The replay coordinator wants to apply remote events to cache without invalidating-and-refetching, to avoid the cross-entity refetch race that would re-introduce banner flicker. This requires the SSE payload to carry post-reanalyze state (the new `analysis_inputs_hash`, updated indices). The payload size grows; under high-traffic conditions this may matter. Measure in the M2 spike; if frame size becomes a concern, fall back to a more compact "you should refetch the following keys" payload and accept the brief banner-gating window.

5. **`useReanalyzeExposure` ergonomics.** The button in [StaleIndicesBanner.tsx:21](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx) currently uses `useMutation().isPending`. Under the queue, this becomes `useMutationState({ filters: { mutationKey, status: 'pending' } }).length > 0`. Verify ergonomics in the M2 implementation; if the boilerplate gets ugly, extract a `useQueueOpStatus` helper hook.

6. **Idempotent-response cache TTL tuning.** The 1-hour TTL is generous. Real retry budgets are seconds-to-minutes. A shorter TTL (5-15 minutes) reduces storage. Decide based on observed retry patterns post-M0.

7. **`MutationCache.getAll()` ordering reliability.** TanStack Query's `MutationCache` stores mutations in a JS `Set`; insertion order is preserved per ECMAScript spec. The spec depends on this. If a future TanStack version changes the underlying storage to something with different iteration semantics, the replay coordinator's queue ordering breaks silently. Mitigation: M1 includes an explicit test that asserts insertion order from `getAll()` after a sequence of `mutate()` calls; this test fails loudly on regression.

## Fallback triggers

Pre-commitments that make the architectural judgment falsifiable in implementation. If any fires, fall back to the alternative described.

- **If any mutator needs to inspect *other* pending mutations to compute its effect.** That is the abstraction-leak signal — the framework's `(op, baseState) → (cacheEffect, request)` shape no longer holds because mutators have hidden dependencies on each other. Fall back to per-mutation optimistic updates via `onMutate` without queue infrastructure. Lose the replay-on-remote-event property; keep autoReanalyze elimination as the primary win.

- **If M0's `analyze_exposure!` fast-skip refactor doesn't reduce the no-change path to microseconds** (target: <100µs total wall-clock; measured via M0 test that runs N curations producing zero net effect and asserts P99 wall-clock). The autoReanalyze elimination story depends on the no-change path being effectively free. If the fast-skip can't deliver this, M2's synchronous-reanalyze pattern relocates latency rather than eliminating it; revert to client-side autoReanalyze and treat the queue as cache-merging-only.

- **If the response-body cache produces user-visible drift** between rehydrated retries and current server state more than a few times per long session (a deployed-then-upgraded-mid-session response shape change is the canonical case). Schema-versioning the cached responses or scoping retries by deploy-version becomes necessary. Diagnostic: M0 includes telemetry on cache-hit retries showing whether the response was generated under the current server version.

- **If two-context Playwright tests in M2 reveal replay-rerun producing user-visible artifacts** (cache flicker, ordering surprises) at rates greater than a small handful per test session, the cache-merged optimistic design (Q7) is exposing seams. Consider whether the M2 failure-class indicators are sufficient to mask, or whether the design needs revision (e.g., a "freeze cache for ~50ms after remote event during pending op" debounce in the replay coordinator).

These triggers are diagnostic, not deadlines. Each names what would falsify a load-bearing assumption; if one fires, the spec returns to design conversation rather than continuing through implementation.

## Out of scope

- Cross-session undo / Cmd+Z (deferred; primitives in place).
- `Last-Event-ID` server-side SSE replay (made unnecessary by HTTP-authoritative confirmation).
- Offline editing / localStorage queue persistence (would require a different design).
- Time-travel UI / state reconstruction at arbitrary historical event id.
- ML-driven conflict resolution.
- Authentication tightening beyond `X-Username` (orthogonal).
- `If-Match` headers + 409 retry on delta-shaped routes (R5b — redirected; replay-as-rerun handles the conflict-resolution problem with different tradeoffs).

## Further reading

- [`docs/event-log.md`](../event-log.md) — `apply_event!` dispatcher contract, hash invariants, SSE semantics. The queue's server-side primitives compose with these unchanged.
- [`docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md`](2026-05-01-multiplayer-instrumentation-design.md) — original Plan 7 design. R0–R5a are foundation for this spec; R5b is redirected.
- [`docs/superpowers/specs/2026-05-02-sse-client-id-design.md`](2026-05-02-sse-client-id-design.md) — PR #32, per-tab `client_id` for multi-tab routing. Direct prerequisite for the queue's SSE confirmation primitive.
- [`packages/HimalayaUI/frontend/src/queries.ts`](../../packages/HimalayaUI/frontend/src/queries.ts) — current mutation surface; the M2 migration target.
- [`packages/HimalayaUI/src/events.jl`](../../packages/HimalayaUI/src/events.jl), [`packages/HimalayaUI/src/actions.jl`](../../packages/HimalayaUI/src/actions.jl) — server-side primitives that M0 extends.
- TanStack Query docs on [optimistic updates](https://tanstack.com/query/v5/docs/framework/react/guides/optimistic-updates) and [`MutationCache`](https://tanstack.com/query/latest/docs/reference/MutationCache) — the canonical primitives the framework layers on.
- [Linear's sync engine talk](https://linear.app/blog/scaling-the-linear-sync-engine) and Matthew Weidner's [survey of sync without CRDTs](https://mattweidner.com/2024/06/04/server-reconciliation.html) — prior art for the pattern this spec instantiates.
