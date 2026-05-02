# Multiplayer + Instrumentation Foundation — Design Spec

**Date:** 2026-05-01
**Status:** Draft

## Context

HimalayaUI was designed for single-user curation. The current data model and pipeline were built to *work*, not to support concurrent editing or to capture user behavior as queryable evidence. Two pressures converge:

1. **Multiplayer.** A small lab team needs to curate the same experiment together — see each other's peak edits, index confirmations, and reject reasons in near real time. Today, two browsers on the same DB can clobber each other's writes silently and never see each other's changes without a manual refresh.

2. **Parameter justification.** Many of Himalaya's analysis parameters are reasonable guesses on small evidence: `MEMBER_REATTACH_RELTOL = 0.05` ([pipeline.jl:32](../../packages/HimalayaUI/src/pipeline.jl)), `SNAP_TOL = 0.0025` ([speculative.jl:11](../../packages/HimalayaUI/src/speculative.jl)), the 0.98 R² hard gate in PhasePanel, the coverage × consistency formula in `Himalaya.score()`, the persistence/sharpness/kneedle blend in `findpeaks`, the recently-tuned high-q trim and prominence floor. Each is a hypothesis the system can't currently retire with evidence because user behavior isn't durably captured in queryable form.

The multiplayer requirement and the instrumentation requirement turn out to share the same root cause: **user intent and machine output are entangled in the same tables**, so neither can be observed cleanly across edits, audited historically, nor reconciled across concurrent actors. The right response is a series of layered refactors that disentangle them, after which both multiplayer and instrumentation become straightforward additions on top of clean foundations.

## Goals

- Multiplayer-safe peak and index curation: two users on the same exposure see each other's edits without manual refresh, with predictable conflict resolution.
- Durable user curation that survives reanalysis without snapshot/restore dances.
- Memoized analysis: `findpeaks` and `indexpeaks` are skipped when their inputs are unchanged.
- A queryable curation event log that captures user intent with enough fidelity to retire parameter guesses with evidence.
- Each refactor delivers value on its own; no step is speculative infrastructure for a later step.

## Non-goals

- Full event-sourcing as the only data model (every read is a fold). We layer materialized views over the log for hot reads.
- Cross-experiment learning, ML-driven phase preference, or active-learning loops. These become *possible* with the log; they are not in scope to build.
- Authentication beyond the existing `X-Username` trust model. Multiplayer makes spoofing more impactful, but it's a separate concern.
- Offline editing / sync-on-reconnect. Clients require server reachability; backfill on reconnect is by re-fetch, not log replay.
- Refactor of `indices`, `index_peaks`, `index_groups`, `index_group_members` into a derived computation (originally proposed as "Refactor 4"). The four-table model isn't broken; replacing it would be a separate decision.

## Pain points in the current design

Three load-bearing examples from `_persist_analysis_inner!` ([pipeline.jl:46-437](../../packages/HimalayaUI/src/pipeline.jl)):

1. **Auto peaks lose IDs on every reanalysis** ([pipeline.jl:113-114](../../packages/HimalayaUI/src/pipeline.jl)). Manual peaks keep IDs. This asymmetry forces:
   - Snapshotting `excluded_qs` by q-value and re-applying after re-detection ([pipeline.jl:108-126](../../packages/HimalayaUI/src/pipeline.jl)).
   - Snapshotting speculative `index_peaks` by q-value and re-resolving ~150 lines later ([pipeline.jl:69-75, 189-371](../../packages/HimalayaUI/src/pipeline.jl)).
   - Custom-group member re-attachment by `(phase, basis)` semantic identity ([pipeline.jl:388-423](../../packages/HimalayaUI/src/pipeline.jl)).

2. **`mark_all_indices_stale!` fires on every peak edit** ([routes_peaks.jl:9-16, 38, 78, 117](../../packages/HimalayaUI/src/routes_peaks.jl)). The flag exists because indices are not eagerly recomputed on the server side. The frontend works around this by chaining peak op → reanalyze → invalidate ([queries.ts:103-114, 116-152](../../packages/HimalayaUI/frontend/src/queries.ts)) — every peak edit is two round trips.

3. **`peaks` table mixes machine output with user curation** ([db.jl:66-75](../../packages/HimalayaUI/src/db.jl)). `source='auto'` rows are ephemeral, `source='manual'` rows are durable, `excluded` overlays user intent on machine rows. Three lifetimes, one table, with the snapshot/restore dance in `_persist_analysis_inner!` to keep them coordinated.

Two separate concurrency hazards:

4. **`ensure_custom_group!` has a TOCTOU race** ([routes_analysis.jl:14-40](../../packages/HimalayaUI/src/routes_analysis.jl)) — SELECT-then-INSERT under concurrent route handlers can produce duplicate custom groups for one exposure. Today's single-user usage almost never triggers this; multiplayer makes it common.

5. **`PATCH /api/exposures/:id/select` is sample-scoped** ([routes_exposures.jl:95-114](../../packages/HimalayaUI/src/routes_exposures.jl)) — clears `selected = 0` across all exposures in the sample, then sets one. Two concurrent selections in the same sample silently clobber.

A capability gap:

6. **`user_actions` is a write-only audit table** ([db.jl:123-131, actions.jl:30-42](../../packages/HimalayaUI/src/db.jl)). It already records `add_peak`, `exclude_peak`, `confirm_index`, `delete_speculative`, `set_status`, and ~15 other actions, but it's never read in the data path. The bones of an event log already exist; they need to be elevated to source-of-truth status with structured payloads.

## Design overview

Five layered refactors, executed in order. Each is a small, independently shippable change.

| # | Refactor | Independent value | Multiplayer enabler? | Instrumentation enabler? |
|---|----------|-------------------|---------------------|-------------------------|
| 0 | Pre-existing race fixes | Closes data-corruption bug | Required prerequisite | — |
| 1 | Diff-based reanalysis preserves auto peak IDs | Simpler pipeline, faster reanalysis, fewer dangling-reference edge cases | — | — |
| 2 | Separate `auto_peaks` from `peak_curations` | Eliminates snapshot/restore dance | Curations become naturally append-friendly | Foundation for log payloads |
| 3 | Content-hash memoization | Skip recompute when inputs unchanged | Hash serves as `If-Match` token | — |
| 4 | Promote `user_actions` to `curation_events` | Audit + replay capability | SSE wire format = log entries | Primary deliverable |
| 5 | SSE multiplayer + 409-on-conflict | The actual multiplayer feature | Primary deliverable | — |

Refactors 0–3 are correctness/clarity improvements that pay for themselves before multiplayer ships. Refactor 4 is the instrumentation deliverable. Refactor 5 is the multiplayer feature, riding on top of foundations laid by 1–4.

---

## R0: Pre-existing race fixes

Before any concurrent-edit infrastructure, fix two bugs that exist today and become more frequent under multiplayer.

### R0.1: Partial unique index on `index_groups`

```sql
CREATE UNIQUE INDEX IF NOT EXISTS idx_one_custom_group_per_exposure
  ON index_groups(exposure_id) WHERE kind = 'custom';
```

Closes the TOCTOU race in `ensure_custom_group!` ([routes_analysis.jl:14-40](../../packages/HimalayaUI/src/routes_analysis.jl)). After this, the function's INSERT either succeeds or fails on the constraint; the SELECT-then-INSERT becomes safe under any number of concurrent handlers.

Add to `migrate_schema!` ([db.jl:142-162](../../packages/HimalayaUI/src/db.jl)) so existing DBs heal on next `open_db`.

### R0.2: `selected` becomes LWW-acceptable

`PATCH /api/exposures/:id/select` ([routes_exposures.jl:95-114](../../packages/HimalayaUI/src/routes_exposures.jl)) is intentionally last-writer-wins under multiplayer: "Bob clicked exposure 7 a millisecond after Alice clicked 6" → 7 is selected, period. The SSE channel (R5) keeps both users' UIs in sync. No code change required, but document the contract: `selected` is sample-scoped client state that does not participate in conflict resolution.

---

## R1: Diff-based reanalysis preserves auto peak IDs

### Current behavior

`_persist_analysis_inner!` ([pipeline.jl:46-437](../../packages/HimalayaUI/src/pipeline.jl)) deletes all `peaks` rows where `source='auto'` and re-inserts fresh. Manual peaks survive. The downstream consequences (~150 lines) work around the resulting ID churn:

- Auto peaks' `excluded` flag is preserved by snapshotting q-values and re-applying ([pipeline.jl:108-126](../../packages/HimalayaUI/src/pipeline.jl)).
- Speculative indices' `index_peaks.peak_id` references would dangle, so they're snapshotted and re-resolved by q-value match ([pipeline.jl:69-75, 189-371](../../packages/HimalayaUI/src/pipeline.jl)).
- Custom-group `index_group_members` would dangle (auto indices also get rebuilt), so they're re-attached by `(phase, basis)` similarity within `MEMBER_REATTACH_RELTOL = 0.05` ([pipeline.jl:32, 388-423](../../packages/HimalayaUI/src/pipeline.jl)).

### Refactored behavior

Replace delete-and-rebuild with diff-and-update:

1. Fetch existing `auto_peaks` for the exposure.
2. Run `findpeaks` on the trace.
3. Match new findpeaks output to existing rows by q within tolerance (`max(EXCLUDE_TOL, abs(q) * 0.001)` — same tolerance as today's snapshot logic at [pipeline.jl:118-126](../../packages/HimalayaUI/src/pipeline.jl)).
4. UPDATE matched rows in place (intensity, prominence, sharpness change; id and `excluded` are preserved).
5. INSERT unmatched new peaks.
6. DELETE auto peaks that no longer correspond to any new detection.

Auto peak IDs survive across reanalyses. `index_peaks` references stay valid for both speculative and auto indices. The speculative-reattachment branch in `_persist_analysis_inner!` reduces from ~180 lines to verifying that each speculative's `index_peaks` rows still point at extant peaks (DELETE the row if not), recomputing basis/r²/score from the surviving members, and falling through to the existing auto-discovery pass for unfilled ratio positions.

For auto **indices**, IDs do not need to be preserved — indices are recomputed from peaks, and their identity isn't meaningful across input changes. The custom-group reattachment-by-(phase, basis) logic stays. (Could be made identity-preserving in a follow-up by matching old auto-index rows to new ones by `(phase, basis)` proximity at update time, but not necessary for R1.)

### Files changed

- [packages/HimalayaUI/src/pipeline.jl](../../packages/HimalayaUI/src/pipeline.jl): rewrite the auto-peak section of `_persist_analysis_inner!` (the DELETE/INSERT block at lines 113-140); the speculative re-resolution block (lines 189-371) shrinks dramatically.

### Test impact

The existing testsets in `test_pipeline.jl` ([packages/HimalayaUI/test/test_pipeline.jl](../../packages/HimalayaUI/test/test_pipeline.jl)) — `"persist_analysis!"`, `"persist_analysis! preserves custom-group members across reanalysis"`, `"analyze_exposure! incorporates manual peaks into candidate indices"`, `"analyze_exposure! ignores excluded auto peaks when scoring candidates"` — must all stay green. Add a new testset asserting auto-peak IDs are preserved across an `analyze_exposure!` round trip when the trace is unchanged.

### Independent of multiplayer

This refactor delivers value on a single-user system: simpler code, faster reanalysis (no peaks DELETE/INSERT round trip), fewer dangling-reference edge cases, and identity-preserving behavior that matches what users intuitively expect.

---

## R2: Separate user curation from machine output

### New schema

Replace the single `peaks` table ([db.jl:66-75](../../packages/HimalayaUI/src/db.jl)) with two tables differentiated by lifetime:

```sql
CREATE TABLE auto_peaks (             -- machine output, regenerated by findpeaks
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL
);

CREATE TABLE peak_curations (         -- user intent, durable
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
    -- For 'exclude': matched against auto_peaks by tol; for 'add': the user's q.
    q           REAL NOT NULL,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- NOTE: sharpness is intentionally NOT stored here. It depends on the trace
-- (`Himalaya.sharpness(I)`); persisting it on the curation row decouples it
-- from the trace, which would make the R3 `analysis_inputs_hash` lie when
-- the trace bytes change but the curation q stays put. Sharpness for `add`
-- curations is re-derived inside `effective_peaks` on every analyze.
```

The `effective_peaks` set fed to `Himalaya.indexpeaks` is computed at analysis time. Sharpness for `add` curations is sampled fresh from `Himalaya.sharpness(I)` on every call, so it stays consistent with the trace:

```julia
function effective_peaks(db, exposure_id, q_grid, I)
    auto     = query("SELECT id, q, sharpness FROM auto_peaks WHERE exposure_id = ?", [exposure_id])
    excludes = query("SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude'", [exposure_id])
    adds     = query("SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add'",  [exposure_id])
    kept = filter(p -> !any(e -> within_tol(e.q, p.q), excludes), auto)
    sharp_full = Himalaya.sharpness(I)  # sampled fresh from current trace
    adds_with_sharp = [(id=a.id, q=a.q, sharpness=sharp_full[argmin(abs.(q_grid .- a.q))]) for a in adds]
    sort_by_q(vcat(kept, adds_with_sharp))
end
```

Index_peaks rows reference either an `auto_peaks.id` or a `peak_curations.id`; a `peak_kind` column on `index_peaks` (added in R2) disambiguates. (Earlier drafts encoded the distinction as a sign on the integer id; that was reverted in favor of the explicit column — arithmetic-on-id is the kind of clever encoding that breaks the moment a downstream query forgets the convention.)

### Pipeline simplification

`analyze_exposure!` ([pipeline.jl:473-533](../../packages/HimalayaUI/src/pipeline.jl)) becomes:

```julia
function analyze_exposure!(db, exposure_id, analysis_dir)
    q, I, σ = load_dat(...)
    new_peaks = Himalaya.findpeaks(q, I, σ)
    diff_update_auto_peaks!(db, exposure_id, new_peaks)  # R1

    eff = effective_peaks(db, exposure_id)
    candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
    group = auto_group(candidates)
    persist_indices_and_groups!(db, exposure_id, candidates, group)
end
```

The snapshot/restore dance disappears: there are no excluded q-values to snapshot (curations are in their own table, untouched by analysis), no manual peaks to merge in (they're already in `peak_curations`), no `index_peaks` re-resolution (peak ids preserved by R1). The current 380-line `_persist_analysis_inner!` shrinks substantially.

### Route handler simplification

`POST /api/exposures/:id/peaks` ([routes_peaks.jl:28-47](../../packages/HimalayaUI/src/routes_peaks.jl)):
- Today: INSERT into `peaks` with `source='manual'`, call `mark_all_indices_stale!`, log action.
- After: INSERT into `peak_curations` with `kind='add'`, log action. (Mark-stale becomes derivable — see R3.)

`PATCH /api/peaks/:id { excluded: true }` ([routes_peaks.jl:49-91](../../packages/HimalayaUI/src/routes_peaks.jl)):
- Today: UPDATE `peaks.excluded`, mark all indices stale, log action.
- After: INSERT into `peak_curations` with `kind='exclude'`, log action. The route's URL shape changes — it's no longer "patch this peak"; it's "create a curation about this peak." Migrate to `POST /api/exposures/:id/curations { kind: 'exclude', target_peak_id: 123 }` for symmetry, or keep the existing URL and let the handler resolve peak → exposure → INSERT curation.

`DELETE /api/peaks/:id` ([routes_peaks.jl:93-123](../../packages/HimalayaUI/src/routes_peaks.jl)):
- Today: DELETE from `peaks`, DELETE from `index_peaks`.
- After: for `auto_peaks` ids, this becomes "INSERT exclude-curation" (auto peaks shouldn't be hard-deleted). For manual entries, it's DELETE the corresponding `peak_curations` row with `kind='add'`.

### Frontend impact

The mutation hooks in [queries.ts:116-152](../../packages/HimalayaUI/frontend/src/queries.ts) call `addPeak`, `removePeak`, `setPeakExcluded`. Their internal API calls change slightly (URL paths if migrating to `/curations`); their hook contracts are unchanged. The peak list returned by `GET /api/exposures/:id/peaks` is the joined `effective_peaks` view — frontend doesn't see the split.

### Migration

```sql
-- 1. Create new tables (above). The disambiguation column on index_peaks is
--    added by table rebuild, NOT by ALTER TABLE — see plan R2.1 step 2's
--    `migrate_r2_widen_index_peaks_pk!`. SQLite cannot widen a PRIMARY KEY
--    via ALTER, and ADD COLUMN ... CHECK doesn't enforce against pre-existing
--    rows; the rebuild handles both. Old rows are copied with peak_kind='auto'
--    before the manual-peak repoint in step 4 flips matching rows to 'curation'.

-- 2. Backfill auto_peaks. peaks PKs were AUTOINCREMENT so ids are stable —
--    existing index_peaks.peak_id values for auto peaks remain valid as-is.
INSERT INTO auto_peaks (id, exposure_id, q, intensity, prominence, sharpness)
SELECT id, exposure_id, q, intensity, prominence, sharpness
FROM peaks WHERE source = 'auto';

-- 3. Backfill exclusion curations.
INSERT INTO peak_curations (exposure_id, kind, q, created_by)
SELECT exposure_id, 'exclude', q, NULL
FROM peaks WHERE source = 'auto' AND excluded = 1;

-- 4. Backfill addition curations and capture old→new id mapping.
--    Important: do this row-by-row in code, not via INSERT-SELECT, so we can
--    capture each new peak_curations.id and use it to repoint index_peaks.
--    For each `peaks` row where source='manual':
--      - INSERT INTO peak_curations (exposure_id, kind, q, created_by) ...
--      - capture new_curation_id = last_insert_rowid()
--      - UPDATE index_peaks SET peak_id = new_curation_id, peak_kind = 'curation'
--          WHERE peak_id = old_peak_id

-- 5. Verify no index_peaks rows reference now-defunct peaks.id values.
--    Any orphans are bugs; fail loudly rather than silently dropping.

-- 6. DROP TABLE peaks.
```

This preserves user-built speculative indices that referenced manual peaks. An earlier draft of the spec deleted the orphaned `index_peaks` rows and counted on reanalysis to re-establish them; that broke speculatives whose anchor or members were manual peaks (the speculative might no longer score well enough to recover under the new effective-peaks shape, silently losing the user's hand-built indexing). The repoint approach makes the migration data-preserving.

`created_by` for backfilled curations is NULL — historical attribution isn't recoverable from `peaks`. Future curations populate it from `X-Username` via `get_or_create_user!` ([actions.jl:16-23](../../packages/HimalayaUI/src/actions.jl)).

### Files changed

- [packages/HimalayaUI/src/db.jl](../../packages/HimalayaUI/src/db.jl): new tables, migration in `migrate_schema!`.
- [packages/HimalayaUI/src/pipeline.jl](../../packages/HimalayaUI/src/pipeline.jl): `analyze_exposure!` simplified; `_persist_analysis_inner!` reduced.
- [packages/HimalayaUI/src/routes_peaks.jl](../../packages/HimalayaUI/src/routes_peaks.jl): all four route handlers rewrite curation operations against `peak_curations`.
- [packages/HimalayaUI/src/json.jl](../../packages/HimalayaUI/src/json.jl): the joined-view shape for `GET /api/exposures/:id/peaks` may need a serialization helper.
- [packages/HimalayaUI/test/test_pipeline.jl](../../packages/HimalayaUI/test/test_pipeline.jl), [test_routes_peaks.jl](../../packages/HimalayaUI/test/test_routes_peaks.jl): rewrite testsets that depend on `peaks(source, excluded)` columns.

### Test impact

Existing testsets covering "incorporates manual peaks" and "ignores excluded auto peaks" must keep green semantics, expressed against the new schema. Add testsets covering: round-trip survival of curations across `analyze_exposure!`, idempotence of double-exclusion (two `exclude` rows for same q resolve to one logical exclusion), and explicit deletion of an `add`-curation row.

---

## R3: Content-hash memoization

### Two new columns on `exposures`

```sql
ALTER TABLE exposures ADD COLUMN trace_hash TEXT;
ALTER TABLE exposures ADD COLUMN analysis_inputs_hash TEXT;
```

- `trace_hash`: SHA-256 over the bytes of the .dat file (computed during `load_dat`).
- `analysis_inputs_hash`: SHA-256 over the canonically-encoded `effective_peaks` set (sorted by q, each peak as `(q::Float64, sharpness::Float64)` reinterpreted to bytes).

### Pipeline becomes hash-guarded

```julia
function analyze_exposure!(db, exposure_id, analysis_dir)
    trace_path = resolve_trace_path(...)
    new_trace_hash = sha256_file(trace_path)

    if new_trace_hash != current_trace_hash(db, exposure_id) || !auto_peaks_exist(db, exposure_id)
        q, I, σ = load_dat(trace_path)
        peaks_result = Himalaya.findpeaks(q, I, σ)
        diff_update_auto_peaks!(db, exposure_id, peaks_result)
        update_trace_hash!(db, exposure_id, new_trace_hash)
    end

    eff = effective_peaks(db, exposure_id)
    new_inputs_hash = hash_peak_set(eff)

    if new_inputs_hash != current_inputs_hash(db, exposure_id) || !indices_exist(db, exposure_id)
        candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
        group = auto_group(candidates)
        persist_indices_and_groups!(db, exposure_id, candidates, group)
        update_inputs_hash!(db, exposure_id, new_inputs_hash)
    end
end
```

### Three properties this delivers at once

**1. Skip recomputation when inputs are unchanged.** Re-analyzing an exposure whose .dat file hasn't changed and whose curation set hasn't changed runs zero analysis code. Today, every reanalysis re-runs `findpeaks` (persistence + sharpness + kneedle) and `indexpeaks` (search across all phases) regardless of whether the result will change.

**2. Idempotent operations don't trigger spurious recomputes.** User excludes peak X then re-includes it → `effective_peaks` is bit-identical to before → hash matches → no recompute, no SSE storm. A counter would have bumped twice; clients would have refetched twice for no reason.

**3. The hash *is* the version token.** `analysis_inputs_hash` serves as the `If-Match` header value for delta-shaped mutations (R5). Two clients converging on the same effective state from different curation paths get the same hash and don't conflict on retry.

### Staleness becomes derivable

The current `mark_all_indices_stale!` mechanism ([routes_peaks.jl:9-16](../../packages/HimalayaUI/src/routes_peaks.jl)) and the `indices.status='stale'` value go away. An index is stale iff the `analysis_inputs_hash` recorded with the index doesn't match the exposure's current `analysis_inputs_hash`. The `StaleIndicesBanner` reads this on render — no UPDATE writes on every peak edit.

Add a column to `indices`:

```sql
ALTER TABLE indices ADD COLUMN inputs_hash TEXT;
```

Set on insert during analysis. Drop the `'stale'` value from `indices.status`'s implicit enum (the column has no CHECK constraint today; just stop writing the value).

### Prior art

`image_version_token` ([image.jl:17-32](../../packages/HimalayaUI/src/image.jl)) already uses the content-version-as-cache-key pattern for detector TIFFs. Combined with `Cache-Control: max-age=31536000, immutable` ([routes_exposures.jl:57-62](../../packages/HimalayaUI/src/routes_exposures.jl)), it lets the browser cache forever while picking up changes when the source file or processing code is bumped. R3 extends the same pattern to analysis state.

### Hash function choice

SHA-256 is overkill for collision resistance at this scale; a 64-bit hash (`Base.hash` over the canonical encoding) would suffice. SHA-256 is preferred for stability across Julia versions (`Base.hash` is documented as not stable across versions), reproducibility under DB inspection, and consistency with file-content hashing. Performance is irrelevant — the trace is hashed once per analyze, the peak set has at most ~30 entries.

### Files changed

- [packages/HimalayaUI/src/db.jl](../../packages/HimalayaUI/src/db.jl): two columns added in migration.
- [packages/HimalayaUI/src/pipeline.jl](../../packages/HimalayaUI/src/pipeline.jl): hash-guard wrapping in `analyze_exposure!`.
- New: a `hash.jl` (or extension to `pipeline.jl`) for `sha256_file` and `hash_peak_set`.
- [packages/HimalayaUI/src/routes_peaks.jl](../../packages/HimalayaUI/src/routes_peaks.jl): `mark_all_indices_stale!` deleted (after the column-driven derivation lands).
- [packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx): condition changes from `index.status === "stale"` to `index.inputs_hash !== exposure.analysis_inputs_hash`.

---

## R4: Promote `user_actions` to `curation_events`

### Why now

`user_actions` ([db.jl:123-131, actions.jl:30-42](../../packages/HimalayaUI/src/db.jl)) already records ~20 distinct curation actions across all routes (`add_peak`, `exclude_peak`, `include_peak`, `remove_peak`, `confirm_index`, `exclude_index`, `create_speculative`, `delete_speculative`, `set_status`, `select_exposure`, `add_tag`, `remove_tag`, `add_message`, `update_sample`, `update_experiment`, `analyze`, `reingest`, `delete_auto_peak_legacy`). It has the right structure and is already wired into every relevant route. It's just write-only — never folded, never queried in the data path.

R4 promotes it to the source of truth for curation state, with `peak_curations` (R2) and `index_group_members` ([db.jl:106-110](../../packages/HimalayaUI/src/db.jl)) becoming materialized views.

### Schema delta

```sql
-- Add structured payload column.
ALTER TABLE user_actions ADD COLUMN payload TEXT;  -- JSON blob

-- Add explicit undo reference.
ALTER TABLE user_actions ADD COLUMN undoes_event_id INTEGER REFERENCES user_actions(id);

-- Index for fold-by-exposure queries.
CREATE INDEX IF NOT EXISTS idx_events_by_exposure
  ON user_actions(entity_type, entity_id, id);
```

Optionally rename `user_actions` to `curation_events`. Pure rename adds clarity but has migration cost; can defer.

### Payload schema discipline

Each event kind carries a payload designed to support both replay and instrumentation. **Designing the payload up front matters more than the storage choice** — fields you forget to record at design time can't be reconstructed later.

Examples (illustrative, not exhaustive):

| Event kind | Payload fields | Supports analysis |
|------------|---------------|-------------------|
| `peak_added` | `q`, `sharpness`, `near_auto_peak_id?`, `predicted_for_phase?`, `predicted_ratio_position?` | "Are users adding peaks at predicted ratio positions of speculative phases?" |
| `peak_excluded` | `auto_peak_id`, `q`, `prominence`, `sharpness` (at exclusion time), `kneedle_position?` | "Which detection-confidence features predict exclusion likelihood?" |
| `peak_unexcluded` | `undoes_event_id`, `peak_curation_id` (the curation row removed) | "How often do users reverse exclusions?" |
| `index_confirmed` | `index_id`, `phase`, `basis`, `score`, `r_squared`, `competing_indices: [{phase, basis, score, r_squared}, ...]` | "Did users prefer the highest-scoring candidate? When they didn't, what did they pick?" |
| `index_unconfirmed` | `index_id`, `undoes_event_id` | "What fraction of confirmations are reverted?" |
| `speculative_created` | `phase`, `ratio_to_peak_id`, `basis`, `score`, `r_squared` | "What manual indexings do users construct? How do their R² values compare to auto candidates?" |
| `speculative_deleted` | `index_id` | "How many speculatives end up rejected?" |
| `analyze_run` | `trace_hash`, `analysis_inputs_hash` (before/after), `findpeaks_skipped: bool`, `indexpeaks_skipped: bool`, `duration_ms` | "How often does memoization save work? How long does analysis take per exposure?" |
| `set_exposure_status` | `prev_status`, `new_status` | "Are exposures usually accepted-then-rejected or vice versa?" |

Critically, `index_confirmed` includes the **competing candidates' scores** at confirmation time, even though that information is currently available in `indices` *for current state*. Indices get rebuilt; the snapshot of "what the user was choosing between" is only durable if captured on the event.

Likewise, `peak_excluded` records prominence/sharpness *at exclusion time*. After the next reanalysis those values may shift; the snapshot preserves what the user was actually looking at when they decided.

### Materialized view contract

**The dispatcher is the only writer to view tables.** Routes do not INSERT/UPDATE/DELETE on `peak_curations`, `index_group_members`, or any view-mat directly. They validate input, then call `apply_event!`, which is the single chokepoint that:

1. Appends the event row to `user_actions`.
2. Dispatches to `update_view_for_event!(db, kind, entity_id, payload, event_id)` to apply the view-side write.
3. Commits the transaction.

This is what makes `rebuild_views_from_log!` a real fold rather than a cosmetic helper: re-running the dispatcher over a log slice produces the same view state that incremental application would.

```julia
function apply_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id=nothing)
    SQLite.transaction(db) do
        event_id = insert_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id)
        update_view_for_event!(db, kind, entity_id, payload, event_id)  # dispatches by kind
        event_id
    end
end
```

`update_view_for_event!` dispatches:
- `peak_added` → INSERT into `peak_curations` with `kind='add'` (route does NOT touch `peak_curations`)
- `peak_excluded` → INSERT into `peak_curations` with `kind='exclude'`
- `peak_unexcluded` → DELETE the corresponding `peak_curations` row (matched by `q` snapshotted from the original `peak_excluded` event referenced via `undoes_event_id`; the route looks up the original event before delegating)
- `index_confirmed` → ensure custom group exists (idempotent under R0.1's unique constraint), INSERT into `index_group_members`
- `index_unconfirmed` → DELETE from `index_group_members`
- `set_exposure_status` → UPDATE `exposures.status` (LWW)
- `analyze_run` → no view update (the analysis pipeline writes `auto_peaks` and `indices` directly; the event is purely instrumentation)

A `rebuild_views_from_log!(db, exposure_id)` function exists for migration, disaster recovery, and replay testing. Tested via property: starting from empty views, applying every event in order via `update_view_for_event!` produces the same state as the live views maintained incrementally by `apply_event!`. **This property test is the contract enforcer** — if a route bypasses `apply_event!` and writes a view directly, the test fails because the log doesn't fold to the same state.

### Indices/groups get the same treatment

`index_groups`, `index_group_members` ([db.jl:97-110](../../packages/HimalayaUI/src/db.jl)) become materialized views of `index_confirmed`/`index_unconfirmed` events. The custom-vs-auto-group distinction stays as an internal detail — `index_confirmed` always means "user pinned this index as their answer," which lives in the custom group; auto groups are still server-derived from `auto_group()`.

The R0.1 partial unique index continues to work; `apply_event!` for `index_confirmed` calls `ensure_custom_group!` inside the same transaction.

### Why the view layer is non-negotiable

Pure event-sourcing (every read is a fold over the log) breaks down at the data path: `analyze_exposure!` reads the effective peak set on every invocation; the speculative reattachment branch reads `index_peaks`; the routes return cached views. Folding events on every read is O(N) per read; for a heavily-curated exposure that's tens to hundreds of events, costing milliseconds per analyze. Acceptable, but less honest than just maintaining the view eagerly. The "log is the truth, view is a cache" model gives the audit/replay benefits without the read-time cost.

### Files changed

- [packages/HimalayaUI/src/db.jl](../../packages/HimalayaUI/src/db.jl): `payload` and `undoes_event_id` columns; index for fold queries.
- [packages/HimalayaUI/src/actions.jl](../../packages/HimalayaUI/src/actions.jl): `log_action!` becomes `apply_event!` (or wraps it). Existing call sites need payload arguments — about 20 routes.
- All `routes_*.jl` files that currently call `log_action!`: extended to pass structured payload.
- [packages/HimalayaUI/src/pipeline.jl](../../packages/HimalayaUI/src/pipeline.jl): emits `analyze_run` events.
- New: `events.jl` with `apply_event!`, `update_view_for_event!`, `rebuild_views_from_log!`.

### Test impact

- Existing route tests covering `log_action!` writes need payload assertions.
- New testset: every event kind round-trips through `rebuild_views_from_log!` to the same view state.
- New testset: undoing an event (`peak_unexcluded` with `undoes_event_id`) restores the view to its pre-undone state.

---

## R5: SSE multiplayer + (optional) 409-on-conflict

With R0–R4 in place, multiplayer is a thin layer on top. **Split into two sub-deliverables** so the cheaper, higher-value half can ship first and the more expensive half can be gated on instrumentation data:

- **R5a — SSE broadcast + LWW everywhere.** Notifies clients of others' edits within ~1s; routes are last-writer-wins on contention. Delivers ~80% of the perceived multiplayer value: "I see Bob's edits without refreshing."
- **R5b — `If-Match` conflict resolution on delta-shaped routes.** Adds optimistic-concurrency 409+retry on the routes where two users editing the same exposure in the same hash window could clobber each other's deltas. Worth doing only if R4's `analyze_run` and curation-event logs show actual conflicts happening at meaningful rates in practice.

This split is consistent with the spec's own "each layer earns its place" principle: instead of building both halves together because they sound like one feature, ship the cheaper one and let the instrumentation answer whether the expensive one is needed. R4 turns "do conflicts happen?" from a guess into a query.

### R5a: Server side (broadcast)

A new module-level event bus, alongside the existing `Ref{Union{SQLite.DB, Nothing}}` pattern in [server.jl](../../packages/HimalayaUI/src/server.jl):

```julia
# Each subscriber holds a stream + a per-subscriber pending-frames queue.
# Per-subscriber queues (rather than one shared queue) avoid the
# "fast subscriber blocks all others" failure mode and let `broadcast_event!`
# enqueue without blocking on the slowest reader.
const SSE_SUBSCRIBERS = Ref{Vector{NamedTuple}}([])  # (stream, pending::Channel{String})
const SSE_LOCK        = ReentrantLock()
const SSE_WAKEUP      = Threads.Condition()         # broadcast notifies; readers wait
```

`broadcast_event!` formats the SSE frame once, then `put!` it onto every subscriber's `pending` channel and `notify` the condition. Each subscriber's reader loop blocks on its own channel + the condition with a 15 s heartbeat timeout (concrete mechanism resolved by the R5.0 spike — see plan).

A new SSE endpoint:

```
GET /api/events
  → text/event-stream
  Headers: Cache-Control: no-cache
           Connection: keep-alive
           X-Accel-Buffering: no       (defeats nginx buffering)
  → Server registers the stream, blocks on SSE_WAKEUP, writes one frame per
    new event. Heartbeat comment frame every 15s so reverse proxies don't
    kill the idle connection.
```

The subscriber loop blocks on `wait(SSE_WAKEUP)` rather than polling — `broadcast_event!` calls `notify(SSE_WAKEUP)` after each commit. (Earlier drafts had a 1s `sleep` poll which floored event latency at one second and pinned a worker thread per subscriber.)

`apply_event!` (R4), after committing its transaction, broadcasts the event to all subscribers:

```
event: curation
data: {"id": 4321, "kind": "peak_added", "entity_type": "exposure",
       "entity_id": 42, "actor": "alice", "ts": "...",
       "payload": {...}, "exposure_analysis_inputs_hash": "abc123..."}
```

The `exposure_analysis_inputs_hash` field is included on events that affect analysis state, letting clients short-circuit invalidation if their cached hash already matches.

**Best-effort delivery.** If the process dies between transaction commit and broadcast, the event is durable in `user_actions` but no SSE frame is sent. Clients reconcile on reconnect via TanStack Query refetch (cheap and correct for steady state). The principled fix — `Last-Event-ID` header with server-side replay from the log — is deferred to a v2 (see Open Questions); not required for R5a to ship.

### R5b: `If-Match` on delta-shaped routes

(Gated on instrumentation data showing actual conflict frequency. May never ship if R4 reveals contention is rare.)


Only the routes that send deltas computed against prior state require `If-Match` checking:

| Route | `If-Match` value |
|-------|-----------------|
| `POST /api/exposures/:id/peaks` (manual peak add) | `analysis_inputs_hash` of the exposure |
| `DELETE /api/peaks/:id` (curation removal of an `add`) | `analysis_inputs_hash` |
| `POST /api/exposures/:id/curations` (R2-shaped exclude/un-exclude) | `analysis_inputs_hash` |
| `POST /api/groups/:id/members` (confirm index) | `analysis_inputs_hash` |
| `DELETE /api/groups/:id/members/:index_id` (unconfirm) | `analysis_inputs_hash` |
| `POST /api/exposures/:id/speculative` | `analysis_inputs_hash` |
| `DELETE /api/indices/:id` (speculative delete) | `analysis_inputs_hash` |
| `POST /api/exposures/:id/analyze` (manual reanalyze) | `analysis_inputs_hash` |

Routes that don't get `If-Match` (LWW-acceptable):
- `PATCH /api/exposures/:id/status` — scalar
- `PATCH /api/exposures/:id/select` — sample-scoped, LWW (R0.2)
- `PATCH /api/samples/:id` — scalar fields
- Any tag add/remove (additive)
- `POST /api/samples/:id/messages` (chat append-only)

On `If-Match` mismatch: HTTP 409 + body `{ current_inputs_hash: "...", current_state: { peaks: [...], curations: [...] } }`.

### Frontend side (R5a)

A new SSE subscriber in `App.tsx`, opens `EventSource('/api/events')`, dispatches each event:

```typescript
useEffect(() => {
  const es = new EventSource('/api/events');
  es.addEventListener('curation', (e) => {
    const event = JSON.parse(e.data);
    // Self-echo filter: skip events authored by this client. NOTE: if two tabs
    // share the same X-Username (same lab user with two browsers open), each
    // tab sees its own edits via the optimistic update path AND filters the
    // server echo. Edits from the *other* tab also get filtered, so the second
    // tab won't update until manual refetch. Acceptable in lab use; if it
    // becomes a problem, switch to a per-tab client id sent on requests.
    if (event.actor === username) return;
    const id = event.entity_id as number;
    queryClient.invalidateQueries({ queryKey: queryKeys.peaks(id) });
    queryClient.invalidateQueries({ queryKey: queryKeys.indices(id) });
    queryClient.invalidateQueries({ queryKey: queryKeys.groups(id) });
    // (the existing invalidateExposure helper at queries.ts:85-94)
  });
  return () => es.close();
}, [username, queryClient]);
```

`useMutationWithRetry` wrapper for delta-shaped mutations:

```typescript
function useMutationWithRetry<TVars, TResult>(opts: {
  mutationFn: (vars: TVars, ifMatch: string) => Promise<TResult>;
  exposureId: number;
}) {
  const qc = useQueryClient();
  return useMutation({
    mutationFn: async (vars: TVars) => {
      const exposure = qc.getQueryData<Exposure>(queryKeys.exposure(opts.exposureId));
      const ifMatch = exposure?.analysis_inputs_hash ?? '';
      try {
        return await opts.mutationFn(vars, ifMatch);
      } catch (e) {
        if (isConflict(e)) {
          await qc.invalidateQueries({ queryKey: queryKeys.exposure(opts.exposureId) });
          await qc.refetchQueries({ queryKey: queryKeys.exposure(opts.exposureId) });
          const fresh = qc.getQueryData<Exposure>(queryKeys.exposure(opts.exposureId));
          return await opts.mutationFn(vars, fresh!.analysis_inputs_hash);
        }
        throw e;
      }
    },
  });
}
```

Per-resource hash retry is attempted at most once. If a second conflict fires, surface to the user as "this exposure changed under you — refresh and try again." This is the rare case where two users contend on the exact same operation simultaneously.

The existing `autoReanalyze` chain ([queries.ts:103-114](../../packages/HimalayaUI/frontend/src/queries.ts)) — currently making every peak edit a two-round-trip operation — is replaced by server-side memoization (R3): the server reanalyzes if the inputs hash changed, skipping the expensive computations otherwise. The client makes one round trip and gets the new state via the SSE broadcast.

### Files changed

- [packages/HimalayaUI/src/server.jl](../../packages/HimalayaUI/src/server.jl): SSE endpoint, subscriber registry, broadcast helper.
- [packages/HimalayaUI/src/events.jl](../../packages/HimalayaUI/src/events.jl) (R4): `apply_event!` calls broadcast after commit.
- All route handlers for delta-shaped routes: extract `If-Match`, compare against `analysis_inputs_hash`, return 409 on mismatch.
- [packages/HimalayaUI/frontend/src/main.tsx](../../packages/HimalayaUI/frontend/src/main.tsx) or [App.tsx](../../packages/HimalayaUI/frontend/src/App.tsx): SSE subscription effect.
- [packages/HimalayaUI/frontend/src/queries.ts](../../packages/HimalayaUI/frontend/src/queries.ts): `useMutationWithRetry` wrapper; remove `autoReanalyze`.
- [packages/HimalayaUI/frontend/src/api.ts](../../packages/HimalayaUI/frontend/src/api.ts): thread `If-Match` header on mutations; thread `analysis_inputs_hash` through query result types.

### Test impact

- Vitest: `useMutationWithRetry` retry path on 409.
- Playwright two-context test: User A and User B both add peaks to the same exposure; both succeed (deltas don't overlap). User A and User B both exclude the same peak: both succeed (idempotent under hashing).
- Julia: `If-Match` mismatch → 409 with current state body.

---

## Consolidated schema changes

```sql
-- R0
CREATE UNIQUE INDEX IF NOT EXISTS idx_one_custom_group_per_exposure
  ON index_groups(exposure_id) WHERE kind = 'custom';

-- R2
CREATE TABLE auto_peaks (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL
);

CREATE TABLE peak_curations (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
    q           REAL NOT NULL,
    -- NOTE: sharpness is intentionally NOT stored here. It depends on the
    -- trace via Himalaya.sharpness(I); persisting it would decouple it from
    -- the trace and make R3's analysis_inputs_hash lie when bytes change but
    -- curation q stays put. effective_peaks samples it fresh at analysis time.
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- index_peaks gains peak_kind to disambiguate auto vs. curation references,
-- AND its PRIMARY KEY widens to (index_id, peak_id, peak_kind). Fresh DBs get
-- both via the CREATE TABLE in SCHEMA. Legacy DBs CANNOT use ALTER TABLE —
-- SQLite doesn't support widening a PRIMARY KEY via ALTER, and ADD COLUMN
-- ... CHECK doesn't enforce the constraint against existing rows. The legacy
-- path is a table rebuild via `migrate_r2_widen_index_peaks_pk!` (plan R2.1
-- step 2), which must run BEFORE `migrate_r2_split_peaks!` because the
-- split's repoint step writes to peak_kind.

-- (peaks table dropped after backfill)

-- R3
ALTER TABLE exposures ADD COLUMN trace_hash            TEXT;
ALTER TABLE exposures ADD COLUMN analysis_inputs_hash  TEXT;
ALTER TABLE indices   ADD COLUMN inputs_hash           TEXT;

-- R4
ALTER TABLE user_actions ADD COLUMN payload          TEXT;  -- JSON
ALTER TABLE user_actions ADD COLUMN undoes_event_id  INTEGER REFERENCES user_actions(id);
CREATE INDEX IF NOT EXISTS idx_events_by_exposure
  ON user_actions(entity_type, entity_id, id);
```

All migrations slot into `migrate_schema!` ([db.jl:142-162](../../packages/HimalayaUI/src/db.jl)) using the existing idempotent ALTER-with-try/catch pattern. The R2 `peaks` → `auto_peaks` + `peak_curations` migration uses the same approach as `migrate_pk_to_autoincrement!` ([db.jl:178-227](../../packages/HimalayaUI/src/db.jl)): create new tables, copy via INSERT-SELECT, drop the old table inside a `SQLite.transaction` block. FK enforcement is toggled via `PRAGMA foreign_keys = OFF` immediately *before* the transaction and restored in a `finally` after — SQLite requires the toggle outside an open transaction.

## Files affected (consolidated)

Backend:
- [packages/HimalayaUI/src/db.jl](../../packages/HimalayaUI/src/db.jl) — schema and migrations
- [packages/HimalayaUI/src/pipeline.jl](../../packages/HimalayaUI/src/pipeline.jl) — diff-based reanalysis, hash-guarded analysis, `effective_peaks` helper
- [packages/HimalayaUI/src/actions.jl](../../packages/HimalayaUI/src/actions.jl) — `log_action!` → `apply_event!` with payload
- [packages/HimalayaUI/src/server.jl](../../packages/HimalayaUI/src/server.jl) — SSE endpoint, broadcast registry
- [packages/HimalayaUI/src/routes_peaks.jl](../../packages/HimalayaUI/src/routes_peaks.jl) — curation-shaped routes, `If-Match` handling, `mark_all_indices_stale!` removed
- [packages/HimalayaUI/src/routes_analysis.jl](../../packages/HimalayaUI/src/routes_analysis.jl) — `ensure_custom_group!` simplified by R0.1, `If-Match` on confirm/unconfirm/speculative routes
- [packages/HimalayaUI/src/routes_exposures.jl](../../packages/HimalayaUI/src/routes_exposures.jl) — `If-Match` on `POST /:id/analyze`; status/select stay LWW
- New: `packages/HimalayaUI/src/events.jl` — `apply_event!`, `update_view_for_event!`, `rebuild_views_from_log!`
- New: hash helpers (could go in `pipeline.jl` or a small `hash.jl`)

Frontend:
- [packages/HimalayaUI/frontend/src/queries.ts](../../packages/HimalayaUI/frontend/src/queries.ts) — `useMutationWithRetry`, removal of `autoReanalyze`, thread `analysis_inputs_hash` on cached entities
- [packages/HimalayaUI/frontend/src/api.ts](../../packages/HimalayaUI/frontend/src/api.ts) — `If-Match` header, hash field in result types
- [packages/HimalayaUI/frontend/src/App.tsx](../../packages/HimalayaUI/frontend/src/App.tsx) or `main.tsx` — SSE subscriber
- [packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx](../../packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx) — staleness derived from hash mismatch

Tests:
- [packages/HimalayaUI/test/test_pipeline.jl](../../packages/HimalayaUI/test/test_pipeline.jl) — testsets adapted to new schema; new ID-preservation testset; new memoization testset
- [packages/HimalayaUI/test/test_routes_peaks.jl](../../packages/HimalayaUI/test/test_routes_peaks.jl) — curation-shaped routes
- [packages/HimalayaUI/test/test_speculative.jl](../../packages/HimalayaUI/test/test_speculative.jl) — speculative reattachment after R1
- New: tests for `apply_event!` view consistency, SSE broadcast, `If-Match` behavior
- [packages/HimalayaUI/frontend/test/](../../packages/HimalayaUI/frontend/test/) — `useMutationWithRetry` unit tests
- [packages/HimalayaUI/frontend/e2e/](../../packages/HimalayaUI/frontend/e2e/) — two-context multiplayer E2E

## Migration plan

Each refactor is its own PR, mergeable independently:

1. **R0** — Two small commits (partial unique index; document `selected` LWW).
2. **R1** — Diff-based reanalysis. Existing pipeline tests stay green; new ID-preservation test added.
3. **R2** — Schema migration + route rewrite. Bigger PR; backfill script tested against fresh and migrated DBs. Migration repoints `index_peaks` rows that referenced manual peaks at the new `peak_curations` ids so user-built speculatives survive.
4. **R3** — Hash columns + memoization. New helpers, no contract changes from outside.
5. **R4** — Promote `user_actions`. Payload backfill is impossible for historical events (data not recorded); new events have full payloads from this PR forward. Routes route view-writes through `apply_event!` exclusively.
6. **R5a** — SSE broadcast + LWW. Spike Oxygen.jl 1.10.x streaming API first (concrete API choice not yet confirmed); commit downstream work only after spike resolves.
7. **R5b** — `If-Match` + 409 retry. Gated on R4 instrumentation data showing actual conflict frequency. May never be needed.

The on-disk DB heals on next `open_db` for any subset of these merged. Partial deployment (e.g., R0+R1+R3 in production while R2 is in review) is safe because each refactor is internally consistent.

## Open questions

1. **Hash function: SHA-256 vs. cheap 64-bit?** SHA-256 is recommended for stability, but `Base.hash` over a canonical encoding would be 100× faster. The peak set hash is already cheap (~30 entries); the trace hash is the only one with measurable cost (~megabyte per file × thousands of exposures across a re-ingest). Proposal: SHA-256 for trace, fast 64-bit for peak set, document the choice.

2. **Event payload completeness review.** Before R4 ships, walk through the parameter-justification questions enumerated in Context and verify each can be answered from the proposed payloads. Adding fields after the fact loses fidelity for historical events.

3. **Hash-as-version-token vs. event-id-as-version-token.** ~~Open.~~ Resolved: hash for delta-mutations (idempotent-equivalence — Alice and Bob converging on the same state both succeed, which matters for retry), event id for SSE cursor (strict serialization, needed for `Last-Event-ID` replay if v2 ships it). Both purposes coexist.

4. **Reconnect backfill.** SSE clients that drop and reconnect will miss events broadcast in the gap. Simplest fix: on reconnect, refetch the affected resources (TanStack Query already supports this). More principled: SSE `Last-Event-ID` header with server-side replay from `user_actions`. Defer the principled version to v2.

5. **Authentication tightening.** `X-Username` is trusted today. Multiplayer doesn't change the trust model but raises the impact of spoofing (Alice's UI shows "Bob excluded peak X" with Bob's name). Out of scope for this spec; flagged for later.

6. **Speculative-index handling under multiplayer.** R1's diff-based reanalysis preserves auto peak ids; speculative `index_peaks` references stay valid. But if Alice's `peak_added` event lands at a predicted ratio position of Bob's existing speculative, Bob's UI suddenly shows the speculative referencing a peak Bob didn't pick. This is correct behavior (the speculative survived reanalysis under new inputs) but warrants a UX affordance. Surface as a notification: "your Pn3m speculative now claims a peak Alice added."

7. **Is the `peaks` → `auto_peaks` + `peak_curations` table split worth the migration cost?** Alternative: keep `peaks` as is, add a `created_by` column for manual peaks, accept the entanglement. Considerations on each side:

   - **R1 alone doesn't fix the snapshot/restore dance.** Even with auto peak ids preserved (R1), `_persist_analysis_inner!` still has to merge manual peaks into the indexpeaks input and re-apply excluded-q snapshots — these live in the pipeline because manual + auto + excluded share one table. Only the table split eliminates the dance.
   - **R4's event payloads reference `peak_curations.id`.** The instrumentation case for R2 is structural: the event log payloads talk about curation rows, which are first-class entities. With the unified `peaks` table, payloads have to reference manual peaks by `peaks.id` and exclusions by `(peak_id, excluded=1)` — workable but harder to fold and harder to query historically.
   - **The `created_by`-only alternative ships faster.** If the only goal were multiplayer correctness, adding `peaks.created_by` and SSE-broadcasting on writes would suffice for R5a.

   Recommend the split. If R2 implementation hits any of the following triggers, fall back to the `peaks.created_by` path instead of pushing through:

   - The R2.2 `_persist_analysis_inner!` rewrite breaks more than 3 of the existing pipeline testsets in ways that don't have an obvious one-line fix.
   - The R2.1 migration backfill produces orphan `index_peaks` rows on a real production-DB sample (not just synthetic fixtures).
   - The cumulative diff for R2.1 + R2.2 exceeds ~600 net lines added (rough proxy for "this turned into a rewrite, not a refactor"). The current pipeline.jl is ~530 lines; a clean R2 should *shrink* it.

   Pre-committing the trigger conditions makes the fallback decision actionable mid-implementation rather than rationalized after the fact. The R1+`created_by` shortcut delivers ~80% of the multiplayer story (R0+R1+R3+R4+R5a all still ship; only R4 payloads become slightly less ergonomic and `_persist_analysis_inner!` keeps its snapshot/restore dance).

## Out of scope

- ML-driven analysis improvements (item 7 from the broader instrumentation discussion). Possible after R4; not committed.
- Cross-experiment learning (priors that span experiments). Possible after R4; not committed.
- Indices as derived computation (eliminate `indices`, `index_peaks` tables). Could be a follow-up after the system has run with R3 memoization for a while.
- Authentication beyond `X-Username`.
- Offline editing / replay-on-reconnect (open question 4 covers a simpler version).
- Time-travel UI (reconstruct an exposure's state at an arbitrary historical event id). The capability is enabled by R4; the UI surfacing is deferred.

## Further reading

- [docs/peak-finding.md](../peak-finding.md) — `findpeaks` design, parameters that R4 instrumentation aims to retire with evidence.
- [docs/scoring.md](../scoring.md) — `score()` formula, weights that R4 will let us validate.
- [docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md](2026-04-22-himalaya-web-app-design.md) — original web app design (single-user assumptions throughout).
- [docs/superpowers/specs/2026-04-23-index-scoring-design.md](2026-04-23-index-scoring-design.md) — prior parameter-tuning decision; an example of the kind of decision R4 evidence would inform.
