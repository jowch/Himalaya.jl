# Multiplayer + Instrumentation Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the five refactors specified in [docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md](../specs/2026-05-01-multiplayer-instrumentation-design.md): close pre-existing concurrency races, simplify the analysis pipeline via diff-based reanalysis and content-hash memoization, separate user curation from machine output, promote `user_actions` to a structured event log, and add SSE-driven multiplayer with optimistic-concurrency conflict resolution.

**Architecture:** Sequential refactors R0 → R5, each independently shippable. Each refactor either changes schema, changes data, or changes code; this plan treats migration as a first-class concern with explicit idempotency, ordering, and rollback constraints.

**Tech stack:** Julia 1.9+, SQLite.jl, Oxygen.jl, HTTP.jl (for SSE), JSON3, existing HimalayaUI conventions; React 18, TanStack Query v5, Vitest, Playwright.

**Read first:**
- [docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md](../specs/2026-05-01-multiplayer-instrumentation-design.md) — the full design rationale
- [CLAUDE.md](../../../CLAUDE.md) — HimalayaUI gotchas (SQLite.jl row materialization, FK enforcement, AUTOINCREMENT contract for mention targets)

---

## Migration Architecture

> **This section governs all DB-touching work in this plan.** Read it before starting any task that touches `db.jl` or any `routes_*.jl` file. Each task's "Migration impact" subsection refers back to the contracts established here.

### Deployment context

- **One central DB** per deployment, resolved by `default_db_path()` in [db.jl:298-303](../../packages/HimalayaUI/src/db.jl) — env `HIMALAYA_DB_PATH` or `~/.himalaya/himalaya.db`. Migrations apply to whatever DB is open at `serve` / `init` / `analyze` startup.
- **`open_db` runs `create_schema!` then `migrate_schema!`** on every connection ([db.jl:315-340](../../packages/HimalayaUI/src/db.jl)). Both must remain idempotent across all DB states this plan introduces.
- **`Manifest.toml` is gitignored** so worktrees re-resolve. New stdlib deps added during this work (e.g., `SHA`) must be added to `Project.toml` `[deps]` and `[compat]` blocks.
- **The DB file is shared across user sessions** (umask 0664 fix in [db.jl:329-338](../../packages/HimalayaUI/src/db.jl)) — migrations may run with multiple readers connected. Each schema change runs inside a single `SQLite.transaction` so partial states are never observable.

### Idempotency contract

Every migration step in this plan must be safe to run on:

1. **Fresh DB** (created by `create_schema!` with current schema)
2. **Pre-this-plan DB** (one or more pre-R0 fields/tables present)
3. **Partially-migrated DB** (some refactors applied, others not)
4. **Already-migrated DB** (all refactors applied; subsequent `open_db` is a no-op)

Concrete patterns:
- `CREATE TABLE IF NOT EXISTS` — idempotent.
- `ALTER TABLE ... ADD COLUMN` — wrap in `try ... catch end` (SQLite errors on duplicate column; the existing pattern at [db.jl:142-160](../../packages/HimalayaUI/src/db.jl) is the precedent).
- `CREATE INDEX IF NOT EXISTS` — idempotent.
- **Table-rebuild migrations** (R2's `peaks` split): use the sentinel pattern from `migrate_pk_to_autoincrement!` ([db.jl:178-227](../../packages/HimalayaUI/src/db.jl)) — detect "is migration needed?" by querying `sqlite_master` or by checking for absence of new tables, return early if not.
- **Data backfills**: each backfill step must be guarded by "have we already run this?" check (typically: target table is empty, source table exists). Re-running on a migrated DB must be a no-op.

### Ordering constraint

Migrations must run in R0 → R1 → R2 → R3 → R4 → R5 order on a per-step basis. The order matters because:
- R2 (peaks split) reads from old `peaks` table, then drops it. R3 hashes the post-R2 effective peak set — running R3's hash code before R2's split would hash the wrong shape.
- R4 (event log) writes payloads that, in R5, are read on SSE broadcast — R5 cannot ship before R4 lands the payload column.
- R5 (`If-Match`) checks `analysis_inputs_hash` from R3 — cannot ship before R3.

`migrate_schema!` should call ordered migration sub-functions:

```julia
function migrate_schema!(db::SQLite.DB)
    migrate_pk_to_autoincrement!(db)        # existing (pre-this-plan)
    migrate_r0_indexes!(db)                 # R0
    migrate_r2_widen_index_peaks_pk!(db)    # R2 — must run before split (repoint
                                             #       step writes peak_kind column)
    migrate_r2_split_peaks!(db)             # R2 (R1 has no schema change)
    migrate_r3_hash_columns!(db)            # R3
    migrate_r4_event_columns!(db)           # R4 (also calls
                                             #     migrate_r4_rebase_entity_type!
                                             #     before creating idx_events_by_exposure)
    # R5 has no schema change — server endpoints + frontend
end
```

Each `migrate_rN_*` is internally idempotent. The order is determined by data dependencies: R2 must complete before R3 because R3 reads `auto_peaks` to compute `analysis_inputs_hash`; R3 must complete before R4 because the event log can reference hashes.

### Rollback story

SQLite ALTER TABLE is **forward-only** in practice — there is no built-in column-drop until SQLite 3.35 and rebuilding tables to remove a column is destructive. The rollback story is therefore:

- **Before any DB-touching task: full DB backup is mandatory.** `cp ~/.himalaya/himalaya.db ~/.himalaya/himalaya.db.pre-rN-backup`.
- **No automated rollback.** If a migration fails mid-deploy, restore from backup and investigate. Do not attempt to manually undo partial state.
- **Each migration commits atomically** (single transaction wrapping the per-step work). A failure inside the transaction leaves the DB unchanged; failures across transactions leave the DB at a consistent intermediate state (some migrations applied, others not), which the next `open_db` will resume from.
- **Migration logging.** Each `migrate_rN_*` function logs `@info` on entry and exit with row counts for backfill steps. Operators can audit `journalctl` (or whatever the daemon log target is) to see what ran.

### DB state matrix (test fixtures)

Every DB-touching task adds at least three test cases (Julia + Vitest where relevant):

| Fixture | Description | What's tested |
|---------|-------------|---------------|
| **Fresh** | `:memory:` DB freshly created via `open_db` | Schema is current; new code paths work |
| **Pre-migration** | DB constructed by recreating the *pre-this-plan* schema, then populated with rows | `migrate_rN_*` correctly transforms historical data |
| **Already-migrated** | A fresh DB; `open_db` called twice | Migration is a no-op the second time (no errors, no row duplication) |

For R2 specifically (table rebuild), add a fourth fixture: a DB with the `peaks` table half-populated mid-transaction is **not** a real failure mode (transactions are atomic), so don't write that test — but do test that aborting a connection mid-`SQLite.transaction` leaves the DB unchanged.

### Deployment coupling between refactors

Most refactors here are independent (each lands behind its own commit cadence), but **R2 and R3 must ship in the same release**. R2's `migrate_r2_split_peaks!` sets `status='stale'` on every existing index *and* deletes the `mark_all_indices_stale!` call sites; R3's hash-mismatch check is what replaces them. Shipping R2 alone leaves `StaleIndicesBanner` permanently visible with no client-visible path to clear it. Treat R2.1 + R2.2 + R3.1 + R3.2 as a single deployment unit. R2.3 (drop legacy `peaks` table) follows separately, after sustained-use confirmation. R0, R1, R4, R5a, R5b can each ship independently of one another (subject to the ordering constraint above).

### Code-side migration

Beyond schema, several refactors change code contracts that consumers depend on:

- **R2:** all consumers of `peaks(source, excluded)` columns need rewriting. `routes_peaks.jl`, `pipeline.jl`, `speculative.jl` — and the frontend's peak list shape (returned by `GET /api/exposures/:id/peaks`) stays the same outwardly: the route handler joins `auto_peaks` and `peak_curations` to produce the legacy shape until R2's frontend pieces ship. This lets backend R2 land without a frontend co-deploy.
- **R3:** `mark_all_indices_stale!` is removed. Routes that called it must drop the call. The `indices.status='stale'` check in `StaleIndicesBanner` becomes a hash-equality check; this is a coordinated frontend+backend change.
- **R4:** `log_action!` becomes `apply_event!` with payload arg. Existing call sites need to either pass a payload or use a thin compat wrapper that defaults `payload=nothing`. Recommended: ship `apply_event!` alongside `log_action!`, migrate call sites in batches, then deprecate `log_action!`.
- **R5:** server adds new endpoints; frontend adds SSE subscriber + `If-Match` retry. Routes accept-but-ignore missing `If-Match` for one release window so old browser tabs don't break, then start enforcing.

### Pre-flight checklist

Before any task that touches `db.jl`:

- [ ] DB backup taken (`cp <db> <db>.backup-<task-id>`)
- [ ] Worktree on the right branch (verify `git status`, `git branch`)
- [ ] `Manifest.toml` copied from main (worktree gotcha — see CLAUDE.md "Himalaya core resolution in worktrees")
- [ ] `Pkg.test("HimalayaUI")` green before touching anything
- [ ] Frontend `npm test` and `npm run build` green

---

## File Map (consolidated)

| Action | Path | Touched in |
|--------|------|------------|
| Modify | `packages/HimalayaUI/src/db.jl` | R0.1, R2, R3, R4 |
| Create | `packages/HimalayaUI/src/events.jl` | R4 |
| Create | `packages/HimalayaUI/src/hash.jl` | R3 |
| Modify | `packages/HimalayaUI/src/pipeline.jl` | R1, R2, R3, R4 |
| Modify | `packages/HimalayaUI/src/actions.jl` | R4 |
| Modify | `packages/HimalayaUI/src/routes_peaks.jl` | R2, R3, R5 |
| Modify | `packages/HimalayaUI/src/routes_analysis.jl` | R0.1, R5 |
| Modify | `packages/HimalayaUI/src/routes_exposures.jl` | R5 |
| Modify | `packages/HimalayaUI/src/server.jl` | R5 |
| Modify | `packages/HimalayaUI/src/speculative.jl` | R1 |
| Modify | `packages/HimalayaUI/src/json.jl` | R2 |
| Modify | `packages/HimalayaUI/Project.toml` | R3 (SHA dep) |
| Modify | `packages/HimalayaUI/test/test_pipeline.jl` | R1, R2, R3 |
| Modify | `packages/HimalayaUI/test/test_routes_peaks.jl` | R2, R5 |
| Modify | `packages/HimalayaUI/test/test_speculative.jl` | R1 |
| Modify | `packages/HimalayaUI/test/test_db.jl` | R0–R4 |
| Create | `packages/HimalayaUI/test/test_events.jl` | R4 |
| Create | `packages/HimalayaUI/test/test_sse.jl` | R5 |
| Modify | `packages/HimalayaUI/frontend/src/queries.ts` | R3, R5 |
| Modify | `packages/HimalayaUI/frontend/src/api.ts` | R3, R5 |
| Modify | `packages/HimalayaUI/frontend/src/App.tsx` | R5 |
| Modify | `packages/HimalayaUI/frontend/src/components/StaleIndicesBanner.tsx` | R3 |
| Modify | `packages/HimalayaUI/frontend/test/*` | R3, R5 |
| Modify | `packages/HimalayaUI/frontend/e2e/*` | R5 |

---

## R0: Pre-existing race fixes

### Task R0.1: Partial unique index on `index_groups`

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

**Migration impact:** Adds one partial unique index. No data migration. Idempotent via `IF NOT EXISTS`. Pre-flight: verify no existing duplicate-custom-group rows would violate the constraint (a SELECT before the CREATE that errors loudly if duplicates exist — those would have to be merged manually before the constraint can land).

- [ ] **Step 1: Write the failing test**

  In `test_db.jl`, add:
  ```julia
  @testset "index_groups partial unique constraint on custom" begin
      mktempdir() do dir
          db = HimalayaUI.open_db(joinpath(dir, "h.db"))
          # set up an exposure
          exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
          s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
          e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

          # First custom group succeeds.
          DBInterface.execute(db,
              "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])
          # Second custom group on same exposure must fail.
          @test_throws SQLite.SQLiteException DBInterface.execute(db,
              "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'custom')", [e_id])

          # Auto groups can multiply (only 'custom' is unique per exposure).
          DBInterface.execute(db,
              "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'auto')", [e_id])
          DBInterface.execute(db,
              "INSERT INTO index_groups (exposure_id, kind) VALUES (?, 'auto')", [e_id])
      end
  end
  ```

- [ ] **Step 2: Run to confirm failure**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  ```

- [ ] **Step 3: Pre-flight check for existing duplicates** (must run before index creation)

  The pre-flight has to fire **before** `CREATE UNIQUE INDEX` runs, otherwise SQLite's own integrity check errors first with a low-context message like "UNIQUE constraint failed: index_groups.exposure_id" — the operator wouldn't know that "merge the duplicate custom groups" is the right next step. Since `create_schema!` is the path that runs the CREATE on fresh and pre-existing DBs alike, the pre-flight goes at the **top of `create_schema!`** (or in `open_db` immediately before `create_schema!` is called):

  ```julia
  # Defensive: if a multiplayer-era duplicate-custom-group row exists in a
  # pre-R0.1 DB, fail loudly with a useful message rather than letting the
  # CREATE UNIQUE INDEX below produce SQLite's terse integrity error.
  # Skipped on truly-fresh DBs (no index_groups table yet).
  has_table = !isempty(Tables.rowtable(DBInterface.execute(db,
      "SELECT 1 FROM sqlite_master WHERE type='table' AND name='index_groups'")))
  if has_table
      dups = Tables.rowtable(DBInterface.execute(db, """
          SELECT exposure_id, COUNT(*) AS n FROM index_groups
          WHERE kind = 'custom' GROUP BY exposure_id HAVING n > 1
      """))
      if !isempty(dups)
          error("DB has duplicate 'custom' index_groups for exposures " *
                join([string(d.exposure_id) for d in dups], ", ") *
                " — manual merge required before idx_one_custom_group_per_exposure can be enforced")
      end
  end
  ```

- [ ] **Step 4: Add the index in `db.jl`**

  Append to `SCHEMA` constant after the `index_groups` table definition:
  ```sql
  CREATE UNIQUE INDEX IF NOT EXISTS idx_one_custom_group_per_exposure
      ON index_groups(exposure_id) WHERE kind = 'custom';
  ```

  No `migrate_schema!` change needed — `IF NOT EXISTS` makes the `create_schema!` path idempotent on existing DBs too. **Verify on a pre-existing DB**: open a real-world DB (or a copy), run `open_db`, confirm the index appears (`SELECT * FROM sqlite_master WHERE name = 'idx_one_custom_group_per_exposure'`).

- [ ] **Step 5: Run tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
  git commit -m "fix(db): enforce one 'custom' index_group per exposure (closes TOCTOU race)"
  ```

### Task R0.2: Document `selected` LWW contract

**Files:**
- Modify: `packages/HimalayaUI/src/routes_exposures.jl` (comment only)
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add comment in route handler**

  In [routes_exposures.jl:95-114](../../packages/HimalayaUI/src/routes_exposures.jl), above the `PATCH /:id/select` handler, add:
  ```julia
  # `selected` is sample-scoped client state — exactly one exposure per sample is
  # marked selected at a time. Under multiplayer this is intentionally LWW: if Alice
  # and Bob click different exposures within the same sample, the later click wins
  # and SSE broadcasts the resulting state. This route does NOT participate in
  # If-Match conflict resolution (see specs/2026-05-01-multiplayer-instrumentation-design.md).
  ```

- [ ] **Step 2: Add to CLAUDE.md gotchas**

  Append to "HimalayaUI gotchas":
  ```markdown
  **`exposures.selected` is sample-scoped LWW.** `PATCH /api/exposures/:id/select` clears `selected = 0` across all exposures in the sample, then sets one. Under multiplayer this is intentional — concurrent selects produce a single resolved value. Don't add If-Match to this route.
  ```

- [ ] **Step 3: Commit**

  ```bash
  git add packages/HimalayaUI/src/routes_exposures.jl CLAUDE.md
  git commit -m "docs: document selected as LWW under multiplayer"
  ```

---

## R1: Diff-based reanalysis preserves auto peak IDs

**Migration impact:** No schema change. The pipeline's behavior changes — auto peak IDs are now stable across reanalyses where they were previously regenerated. **Existing DBs**: on first run after deploy, the next `analyze_exposure!` will go through the new diff path and produce a different distribution of `peaks.id` values than a delete-rebuild would have. This is fine because the IDs were already volatile; nothing external pins them.

The downstream consequence: speculative `index_peaks` references (which currently dangle through reanalysis and are reconstructed by q-value match) now stay valid. The reattachment branch of `_persist_analysis_inner!` simplifies dramatically. **Test the simplification by green tests, not by removing the speculative-recovery branch wholesale** — the recovery code handles cases (peak gone entirely, peak shifted beyond tolerance) that diff-based update doesn't fix.

### Task R1.1: Diff-update helper for auto peaks

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/test/test_pipeline.jl`

- [ ] **Step 1: Write a regression test for ID preservation**

  Add to `test_pipeline.jl`:
  ```julia
  @testset "analyze_exposure! preserves auto peak IDs across reruns" begin
      mktempdir() do dir
          # boilerplate: minimal experiment with one exposure; copy a fixture .dat
          # into analysis_dir; run analyze_exposure! once, capture peak ids.
          db = setup_minimal_experiment(dir)
          exposure_id = first_exposure_id(db)
          HimalayaUI.analyze_exposure!(db, exposure_id, analysis_dir(db, exposure_id))
          ids_before = sort([Int(r.id) for r in HimalayaUI.get_peaks_for_exposure(db, exposure_id)
                                                  if String(r.source) == "auto"])
          # Re-run with no input change.
          HimalayaUI.analyze_exposure!(db, exposure_id, analysis_dir(db, exposure_id))
          ids_after = sort([Int(r.id) for r in HimalayaUI.get_peaks_for_exposure(db, exposure_id)
                                                 if String(r.source) == "auto"])
          @test ids_before == ids_after
      end
  end
  ```

  (You'll need a small fixture-helper module; if one doesn't exist in the test directory yet, factor one out from existing `analyze_exposure!` testsets that already build minimal experiments — see [test_pipeline.jl:146-167](../../packages/HimalayaUI/test/test_pipeline.jl) for the pattern. The `analysis_dir(db, exposure_id)` helper used above is also part of the fixture module — it joins `experiments.analysis_dir` for the exposure's experiment, mirroring the route handlers' `JOIN exposures … samples … experiments` query at [routes_exposures.jl:144-150](../../packages/HimalayaUI/src/routes_exposures.jl). Define it once in the fixture module rather than inlining it per testset.)

- [ ] **Step 2: Run to confirm failure**

  Expected: ids_before ≠ ids_after because today's pipeline deletes and re-inserts.

- [ ] **Step 3: Implement `diff_update_auto_peaks!`**

  In `pipeline.jl`, add a new helper:
  ```julia
  """
  Match new findpeaks output against existing auto rows by q-value within
  tolerance, UPDATE matched rows in place, INSERT unmatched, DELETE orphans.
  Returns Dict{Float64,Int} mapping each new peak's q to its (preserved or new) id.
  """
  function diff_update_auto_peaks!(db::SQLite.DB, exposure_id::Int,
                                    peaks_result::NamedTuple,
                                    I_full::Vector{Float64})
      EXCLUDE_TOL = 1e-6  # matches existing tolerance at pipeline.jl:118
      tol(q) = max(EXCLUDE_TOL, abs(q) * 0.001)

      existing = Tables.rowtable(DBInterface.execute(db,
          "SELECT id, q FROM peaks WHERE exposure_id = ? AND source = 'auto'", [exposure_id]))
      remaining = Set{Int}(Int(r.id) for r in existing)

      q_to_id = Dict{Float64, Int}()
      for i in eachindex(peaks_result.q)
          qval     = peaks_result.q[i]
          full_idx = peaks_result.indices[i]
          intensity, prom, sharp = I_full[full_idx],
                                    peaks_result.prominence[i],
                                    peaks_result.sharpness[i]

          # Find closest existing auto peak within tolerance.
          best_id, best_d = 0, Inf
          for r in existing
              Int(r.id) in remaining || continue
              d = abs(Float64(r.q) - qval)
              if d < best_d && d <= tol(qval)
                  best_d, best_id = d, Int(r.id)
              end
          end

          if best_id != 0
              DBInterface.execute(db,
                  """UPDATE peaks SET q = ?, intensity = ?, prominence = ?, sharpness = ?
                     WHERE id = ?""",
                  [qval, intensity, prom, sharp, best_id])
              delete!(remaining, best_id)
              q_to_id[qval] = best_id
          else
              res = DBInterface.execute(db,
                  """INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source, excluded)
                     VALUES (?, ?, ?, ?, ?, 'auto', 0)""",
                  [exposure_id, qval, intensity, prom, sharp])
              q_to_id[qval] = Int(DBInterface.lastrowid(res))
          end
      end

      # Drop auto peaks that no longer correspond to any new detection. CASCADE
      # the index_peaks rows that referenced them — they were stale anyway.
      for orphan_id in remaining
          DBInterface.execute(db, "DELETE FROM index_peaks WHERE peak_id = ?", [orphan_id])
          DBInterface.execute(db, "DELETE FROM peaks WHERE id = ?", [orphan_id])
      end

      q_to_id
  end
  ```

- [ ] **Step 4: Replace the auto-peak DELETE+INSERT block in `_persist_analysis_inner!`**

  In `pipeline.jl`, locate the existing block at [pipeline.jl:108-140](../../packages/HimalayaUI/src/pipeline.jl) (the snapshot of `excluded_qs`, the `DELETE FROM peaks WHERE source = 'auto'`, and the loop that re-inserts each auto peak with `was_excluded` re-applied). Replace with:
  ```julia
  q_to_peak_id = diff_update_auto_peaks!(db, exposure_id, peaks_result, I_full)
  ```

  The `excluded` flag is preserved by row identity now — no snapshot, no re-application. Manual peak ID merging at [pipeline.jl:153-158](../../packages/HimalayaUI/src/pipeline.jl) stays unchanged: manual peaks were already stable.

- [ ] **Step 5: Verify the speculative re-resolution branch still works**

  Branch [pipeline.jl:189-371](../../packages/HimalayaUI/src/pipeline.jl) snapshots speculative `index_peaks` rows up front and rebuilds them by q-value match. With auto peak IDs now stable, most speculative references are still valid post-update. The recovery code is no longer load-bearing for the *common* case — but it still handles peaks that genuinely went away (auto detection no longer finds them) or shifted beyond tolerance. **Don't delete this branch.** Run the existing testsets that exercise it:
  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI", test_args=["test_speculative.jl", "test_pipeline.jl"])'
  ```

  The "persist_analysis! preserves custom-group members across reanalysis" testset at [test_pipeline.jl:82](../../packages/HimalayaUI/test/test_pipeline.jl) and the speculative-recovery testsets in `test_speculative.jl` must remain green.

- [ ] **Step 6: Run full HimalayaUI test suite**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  ```

  All testsets including the new ID-preservation one must pass.

- [ ] **Step 7: Commit**

  ```bash
  git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_pipeline.jl
  git commit -m "refactor(pipeline): diff-update auto peaks instead of delete-rebuild

  Preserves auto peak IDs across reanalyses, eliminating the q-value snapshot/
  restore dance for excluded flags. Speculative index_peaks references stay
  valid for the common case; the existing recovery branch handles the residual
  cases where peaks genuinely disappeared or shifted beyond tolerance."
  ```

---

## R2: Separate user curation from machine output

**Migration impact:** This is the largest single migration in this plan. New tables `auto_peaks` and `peak_curations`; old `peaks` table is dropped. `index_peaks.peak_id` references must be reconciled. Indices are *not* migrated in place — they're cleared so the next analysis run rebuilds them under the new effective-peaks shape.

**Backup strongly recommended before deploying this task.**

### Task R2.1: Schema additions (without dropping old `peaks`)

Land the new tables alongside the old one. The pipeline still reads/writes `peaks`; this task just makes the new tables exist and adds backfill logic that runs on `open_db`. The actual cutover happens in R2.2.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

- [ ] **Step 1: Add new tables to `SCHEMA`**

  In `db.jl`, append to `SCHEMA` constant:
  ```sql
  CREATE TABLE IF NOT EXISTS auto_peaks (
      id              INTEGER PRIMARY KEY AUTOINCREMENT,
      exposure_id     INTEGER REFERENCES exposures(id),
      q               REAL NOT NULL,
      intensity       REAL,
      prominence      REAL,
      sharpness       REAL,
      findpeaks_index INTEGER  -- grid index from Himalaya.findpeaks; lets
                               -- synthesize_peaks_result reconstruct the exact
                               -- (q, indices, prominence, sharpness) NamedTuple
                               -- without lossy argmin(abs.(q .- pk.q)) lookup.
  );

  CREATE TABLE IF NOT EXISTS peak_curations (
      id          INTEGER PRIMARY KEY AUTOINCREMENT,
      exposure_id INTEGER REFERENCES exposures(id),
      kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
      q           REAL NOT NULL,
      created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
      created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
      -- NOTE: sharpness is intentionally NOT stored. It depends on the trace
      -- via Himalaya.sharpness(I); persisting it here would decouple it from
      -- the trace and make R3's analysis_inputs_hash lie when the trace bytes
      -- change but the curation q stays put. effective_peaks samples it fresh.
  );

  CREATE INDEX IF NOT EXISTS idx_auto_peaks_exposure
      ON auto_peaks(exposure_id);
  CREATE INDEX IF NOT EXISTS idx_peak_curations_exposure
      ON peak_curations(exposure_id);
  ```

- [ ] **Step 2: Add the `index_peaks.peak_kind` column to `SCHEMA`**

  *Edit the existing `CREATE TABLE index_peaks (...)` statement in `SCHEMA` in place* — do not append a second CREATE block (an executing agent who appends would hit a duplicate-statement error on fresh DBs). The post-edit definition:

  ```sql
  -- index_peaks rows can now reference auto_peaks OR peak_curations; the kind
  -- column disambiguates. Existing rows are all 'auto' (manual-peak refs
  -- get repointed during migration).
  CREATE TABLE IF NOT EXISTS index_peaks (
      index_id       INTEGER REFERENCES indices(id),
      peak_id        INTEGER NOT NULL,
      peak_kind      TEXT NOT NULL DEFAULT 'auto'
                     CHECK (peak_kind IN ('auto', 'curation')),
      ratio_position INTEGER,
      residual       REAL,
      PRIMARY KEY (index_id, peak_id, peak_kind)
  );
  ```

  Legacy DBs get the column via ALTER, and the table is rebuilt to widen the PRIMARY KEY (SQLite does not support widening a PK via ALTER — only a rebuild does). The rebuild also picks up the CHECK constraint, which an `ADD COLUMN ... CHECK` does not enforce against pre-existing rows. Append to `migrate_schema!`:

  ```julia
  function migrate_r2_widen_index_peaks_pk!(db::SQLite.DB)
      # Sentinel: skip if peak_kind already in PK (rebuilt previously).
      info = Tables.rowtable(DBInterface.execute(db,
          "SELECT name, pk FROM pragma_table_info('index_peaks')"))
      already_widened = any(r -> String(r.name) == "peak_kind" && Int(r.pk) > 0, info)
      already_widened && return
      # Skip on fresh DBs that already have the new shape via CREATE.
      isempty(info) && return

      DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
      try
          SQLite.transaction(db) do
              DBInterface.execute(db, "ALTER TABLE index_peaks RENAME TO _index_peaks_old")
              # Re-create with the new shape (matches SCHEMA above).
              DBInterface.execute(db, """
                  CREATE TABLE index_peaks (
                      index_id       INTEGER REFERENCES indices(id),
                      peak_id        INTEGER NOT NULL,
                      peak_kind      TEXT NOT NULL DEFAULT 'auto'
                                     CHECK (peak_kind IN ('auto', 'curation')),
                      ratio_position INTEGER,
                      residual       REAL,
                      PRIMARY KEY (index_id, peak_id, peak_kind)
                  )
              """)
              # Old rows are all 'auto' (this runs before migrate_r2_split_peaks!,
              # which is responsible for repointing manual-peak refs).
              DBInterface.execute(db, """
                  INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
                  SELECT index_id, peak_id, 'auto', ratio_position, residual
                  FROM _index_peaks_old
              """)
              DBInterface.execute(db, "DROP TABLE _index_peaks_old")
          end
      finally
          DBInterface.execute(db, "PRAGMA foreign_keys = ON")
      end
  end
  ```

  Wire-in order in `migrate_schema!` matters: this rebuild must run **before** `migrate_r2_split_peaks!`, because the split-peaks migration's repoint step (UPDATE index_peaks SET peak_kind = 'curation') requires the column to exist with the wider PK. Append in this order:
  ```julia
  migrate_r2_widen_index_peaks_pk!(db)  # rebuild with widened PK first
  migrate_r2_split_peaks!(db)            # then repoint manual-peak refs
  ```

  The simpler `ALTER TABLE … ADD COLUMN peak_kind …` form is no longer needed; the rebuild covers both fresh-with-no-old-rows and legacy-with-old-rows cases. Drop any prior ALTER from earlier drafts.

- [ ] **Step 3: Add backfill function with id-preserving repoint**

  The migration must preserve user-built speculative indices that reference manual peaks. Strategy: INSERT each manual peak as a `peak_curations` row, capture `(old_peak_id → new_curation_id)`, then UPDATE matching `index_peaks` rows to point at the new id with `peak_kind='curation'`.

  ```julia
  """
      migrate_r2_split_peaks!(db)

  Backfill `auto_peaks` and `peak_curations` from the legacy `peaks` table,
  repointing `index_peaks.peak_id` for manual-peak references so user-built
  speculatives survive the migration. Idempotent: returns early if `peaks`
  no longer exists.
  """
  function migrate_r2_split_peaks!(db::SQLite.DB)
      # Sentinel: if `peaks` table is gone, migration already ran.
      peaks_exists = !isempty(Tables.rowtable(DBInterface.execute(db,
          "SELECT 1 FROM sqlite_master WHERE type='table' AND name='peaks'")))
      peaks_exists || return

      # Sentinel: if auto_peaks already has rows, we're partway through —
      # the only safe action is to bail and require operator intervention.
      auto_count = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT COUNT(*) AS n FROM auto_peaks"))).n
      if Int(auto_count) > 0
          error("migrate_r2_split_peaks!: auto_peaks already has $auto_count rows " *
                "but peaks table still exists — operator intervention required " *
                "(restore from backup or manually reconcile)")
      end

      DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
      try
          SQLite.transaction(db) do
              # 1. Auto peaks: id preserved (peaks PK was AUTOINCREMENT).
              # findpeaks_index left NULL for legacy rows — synthesize_peaks_result
              # falls back to argmin lookup when NULL; the next analyze run that
              # invokes diff_update_auto_peaks! will populate it.
              DBInterface.execute(db, """
                  INSERT INTO auto_peaks (id, exposure_id, q, intensity, prominence, sharpness, findpeaks_index)
                  SELECT id, exposure_id, q, intensity, prominence, sharpness, NULL
                  FROM peaks WHERE source = 'auto'
              """)

              # 2. Exclusion curations: q-value is the binding key.
              DBInterface.execute(db, """
                  INSERT INTO peak_curations (exposure_id, kind, q, created_by)
                  SELECT exposure_id, 'exclude', q, NULL
                  FROM peaks WHERE source = 'auto' AND excluded = 1
              """)

              # 3. Addition curations: row-by-row to capture old→new id mapping.
              manual_rows = Tables.rowtable(DBInterface.execute(db,
                  "SELECT id, exposure_id, q FROM peaks WHERE source = 'manual'"))
              old_to_new = Dict{Int, Int}()
              for r in manual_rows
                  res = DBInterface.execute(db,
                      """INSERT INTO peak_curations (exposure_id, kind, q, created_by)
                         VALUES (?, 'add', ?, NULL)""",
                      [Int(r.exposure_id), Float64(r.q)])
                  new_id = Int(DBInterface.lastrowid(res))
                  old_to_new[Int(r.id)] = new_id
              end

              # 4. Repoint index_peaks rows whose peak_id was a manual peak.
              for (old_id, new_id) in old_to_new
                  DBInterface.execute(db,
                      """UPDATE index_peaks SET peak_id = ?, peak_kind = 'curation'
                         WHERE peak_id = ?""",
                      [new_id, old_id])
              end

              # 5. Verify no orphan index_peaks rows remain (auto refs survive,
              #    manual refs were just repointed). Any remaining row whose
              #    peak_id doesn't resolve in either table is a bug.
              orphans = Tables.rowtable(DBInterface.execute(db, """
                  SELECT ip.peak_id, ip.peak_kind, ip.index_id
                  FROM index_peaks ip
                  WHERE (ip.peak_kind = 'auto'     AND ip.peak_id NOT IN (SELECT id FROM auto_peaks))
                     OR (ip.peak_kind = 'curation' AND ip.peak_id NOT IN (SELECT id FROM peak_curations))
              """))
              if !isempty(orphans)
                  error("migrate_r2_split_peaks!: $(length(orphans)) orphaned index_peaks " *
                        "rows after repoint — operator intervention required")
              end

              # 6. Mark all indices stale so the next analyze recomputes basis/score
              #    under the new effective_peaks model. R3's hash check supersedes
              #    this once it lands; the explicit flag is a transitional aid.
              #
              # DEPLOYMENT NOTE: ship R2 and R3 together (single release / one
              # `migrate_schema!` invocation). Shipping R2 alone leaves every index
              # tagged 'stale' with no client-visible path to clear it — the
              # frontend's StaleIndicesBanner stays up forever, because R2 deletes
              # `mark_all_indices_stale!` call sites but doesn't yet have R3's
              # hash-based derived-staleness to take their place. R3.2 step 3
              # normalizes 'stale' → 'candidate' and switches the banner to hash
              # comparison, closing the loop.
              DBInterface.execute(db, "UPDATE indices SET status = 'stale'")

              # 7. Drop the old peaks table — fully decomposed and repointed.
              DBInterface.execute(db, "DROP TABLE peaks")
          end
          @info "migrate_r2_split_peaks! complete"
      finally
          DBInterface.execute(db, "PRAGMA foreign_keys = ON")
      end
  end
  ```

- [ ] **Step 4: Wire into `migrate_schema!`**

  Append to `migrate_schema!`:
  ```julia
  migrate_r2_split_peaks!(db)
  ```

  This runs on every `open_db`, but the sentinel makes it a no-op once it's been done.

- [ ] **Step 5: Test on the three required fixtures**

  Add to `test_db.jl`:
  ```julia
  @testset "migrate_r2_split_peaks! on fresh DB is no-op" begin
      mktempdir() do dir
          db = HimalayaUI.open_db(joinpath(dir, "h.db"))
          # Re-run; nothing should change.
          HimalayaUI.migrate_r2_split_peaks!(db)
          @test isempty(Tables.rowtable(DBInterface.execute(db,
              "SELECT 1 FROM sqlite_master WHERE name = 'peaks'")))
      end
  end

  @testset "migrate_r2_split_peaks! on legacy DB backfills" begin
      mktempdir() do dir
          db_path = joinpath(dir, "h.db")
          db = HimalayaUI.open_db(db_path)
          # Recreate the legacy `peaks` table by hand.
          DBInterface.execute(db, "DROP TABLE auto_peaks")
          DBInterface.execute(db, "DROP TABLE peak_curations")
          DBInterface.execute(db, """
              CREATE TABLE peaks (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  exposure_id INTEGER, q REAL, intensity REAL, prominence REAL,
                  sharpness REAL, source TEXT DEFAULT 'auto', excluded INTEGER DEFAULT 0
              )
          """)
          # Set up an exposure
          exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
          s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
          e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
          # Three peaks: one auto kept, one auto excluded, one manual.
          DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.10, 'auto', 0)", [e_id])
          DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.15, 'auto', 1)", [e_id])
          DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q, source, excluded) VALUES (?, 0.20, 'manual', 0)", [e_id])

          # Re-create the (now-missing) destination tables and run migration.
          HimalayaUI.create_schema!(db)
          HimalayaUI.migrate_r2_split_peaks!(db)

          autos = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM auto_peaks"))
          @test length(autos) == 2  # both auto peaks (excluded ones still in auto_peaks)
          curs  = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM peak_curations"))
          @test length(curs) == 2  # one exclude, one add
          @test Set([String(c.kind) for c in curs]) == Set(["exclude", "add"])
          @test isempty(Tables.rowtable(DBInterface.execute(db,
              "SELECT 1 FROM sqlite_master WHERE name = 'peaks'")))
      end
  end

  @testset "migrate_r2_split_peaks! is idempotent" begin
      mktempdir() do dir
          db = HimalayaUI.open_db(joinpath(dir, "h.db"))
          # Open again — should not error or duplicate rows.
          db2 = HimalayaUI.open_db(joinpath(dir, "h.db"))
          @test true  # if no exception thrown, we're good
      end
  end

  @testset "migrate_r2_split_peaks! preserves index_peaks for speculatives anchored on manual peaks" begin
      # The bug this guards against: a user-built speculative index that
      # anchored on a manual peak loses its anchor when migration drops the
      # manual peak. The repoint step in the migration prevents this.
      mktempdir() do dir
          db_path = joinpath(dir, "h.db")
          db = HimalayaUI.open_db(db_path)
          # Construct legacy schema by hand and seed with a speculative
          # referencing a manual peak.
          DBInterface.execute(db, "DROP TABLE auto_peaks")
          DBInterface.execute(db, "DROP TABLE peak_curations")
          DBInterface.execute(db, """
              CREATE TABLE peaks (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  exposure_id INTEGER, q REAL, intensity REAL, prominence REAL,
                  sharpness REAL, source TEXT DEFAULT 'auto', excluded INTEGER DEFAULT 0
              )
          """)
          exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
          s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
          e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
          # Manual peak that the speculative will anchor on.
          res = DBInterface.execute(db,
              "INSERT INTO peaks (exposure_id, q, source) VALUES (?, 0.20, 'manual')", [e_id])
          manual_peak_id = Int(DBInterface.lastrowid(res))
          # Pseudo-speculative index referencing the manual peak.
          res = DBInterface.execute(db,
              "INSERT INTO indices (exposure_id, phase, basis, kind) VALUES (?, 'Pn3m', 0.20, 'speculative')",
              [e_id])
          ix_id = Int(DBInterface.lastrowid(res))
          DBInterface.execute(db,
              "INSERT INTO index_peaks (index_id, peak_id, ratio_position, residual) VALUES (?, ?, 1, 0.0)",
              [ix_id, manual_peak_id])

          # Run migration.
          HimalayaUI.create_schema!(db)
          HimalayaUI.migrate_r2_split_peaks!(db)

          # Assert the speculative still has its anchor — repointed at the
          # new peak_curations.id, not dangling or deleted.
          ip = Tables.rowtable(DBInterface.execute(db,
              "SELECT * FROM index_peaks WHERE index_id = ?", [ix_id]))
          @test length(ip) == 1
          @test String(ip[1].peak_kind) == "curation"
          # The new peak_id should be a valid peak_curations row of the right shape.
          curations = Tables.rowtable(DBInterface.execute(db,
              "SELECT id, kind, q FROM peak_curations WHERE id = ?", [Int(ip[1].peak_id)]))
          @test length(curations) == 1
          @test String(curations[1].kind) == "add"
          @test Float64(curations[1].q) ≈ 0.20
      end
  end
  ```

- [ ] **Step 6: Run tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_db.jl
  git commit -m "feat(db): add auto_peaks + peak_curations tables and migration

  Schema-only change: tables exist alongside legacy 'peaks' until
  the pipeline cutover in the next commit. Backfill is idempotent
  via sentinel checks; aborts on partially-migrated DBs requiring
  operator intervention."
  ```

### Task R2.2: Pipeline cutover

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/src/routes_peaks.jl`
- Modify: `packages/HimalayaUI/src/json.jl`
- Modify: `packages/HimalayaUI/test/test_pipeline.jl`
- Modify: `packages/HimalayaUI/test/test_routes_peaks.jl`

**Migration impact:** This is the code-side cutover. The new `auto_peaks` + `peak_curations` tables are now **read** by the pipeline and routes. The frontend's `GET /api/exposures/:id/peaks` shape stays unchanged: the route handler joins the two tables to produce the legacy peak-row shape (with synthesized `source` and `excluded` fields) so frontend co-deploy is not required.

- [ ] **Step 1: Add `effective_peaks` helper in `pipeline.jl`**

  Sharpness for `add` curations is sampled fresh from the trace via `Himalaya.sharpness(I)` rather than persisted on the curation row. This keeps the R3 `analysis_inputs_hash` honest under trace re-uploads (the hash is computed from the same NamedTuple this returns; if sharpness moved with the trace, hash equality would be a lie).

  Auto vs. curation distinction is carried as an explicit `peak_kind::Vector{Symbol}` field — *not* as a sign on the integer id. Earlier draft had a negative-id sentinel; that was reverted because arithmetic-on-id is exactly the kind of clever encoding that breaks downstream readers.

  ```julia
  """
      effective_peaks(db, exposure_id, q_grid, I) -> NamedTuple

  Compute the effective peak set for analysis. Returns
  `(q::Vector{Float64}, sharpness::Vector{Float64},
    peak_id::Vector{Int}, peak_kind::Vector{Symbol})` sorted by q.
  Auto peaks whose q matches an `exclude` curation are dropped;
  `add` curations are unioned in with sharpness sampled from the
  current trace.
  """
  function effective_peaks(db::SQLite.DB, exposure_id::Int,
                            q_grid::Vector{Float64}, I::Vector{Float64})
      auto = Tables.rowtable(DBInterface.execute(db,
          "SELECT id, q, sharpness FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))
      excludes = Tables.rowtable(DBInterface.execute(db,
          "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude'",
          [exposure_id]))
      adds = Tables.rowtable(DBInterface.execute(db,
          "SELECT id, q FROM peak_curations WHERE exposure_id = ? AND kind = 'add'",
          [exposure_id]))

      tol(q) = max(1e-6, abs(q) * 0.001)
      is_excluded(qv) = any(e -> abs(Float64(e.q) - qv) <= tol(qv), excludes)

      qs        = Float64[]
      shs       = Float64[]
      peak_id   = Int[]
      peak_kind = Symbol[]
      for r in auto
          qv = Float64(r.q)
          is_excluded(qv) && continue
          push!(qs, qv)
          push!(shs, ismissing(r.sharpness) ? 0.0 : Float64(r.sharpness))
          push!(peak_id,   Int(r.id))
          push!(peak_kind, :auto)
      end
      sharp_full = isempty(adds) ? Float64[] : Himalaya.sharpness(I)
      for r in adds
          qv = Float64(r.q)
          push!(qs, qv)
          push!(shs, sharp_full[argmin(abs.(q_grid .- qv))])
          push!(peak_id,   Int(r.id))
          push!(peak_kind, :curation)
      end

      perm = sortperm(qs)
      (q = qs[perm], sharpness = shs[perm],
       peak_id = peak_id[perm], peak_kind = peak_kind[perm])
  end
  ```

- [ ] **Step 2: Rewrite `analyze_exposure!`** (and update `diff_update_auto_peaks!` to target the new tables)

  `diff_update_auto_peaks!` (introduced in R1 against the legacy `peaks` table) is rewritten to (a) read/write `auto_peaks` instead of `peaks(source='auto')`, and (b) populate `findpeaks_index` from `peaks_result.indices[i]` on every INSERT/UPDATE — this is the column R3.2's `synthesize_peaks_result` uses to skip an exact-vs-argmin reconstruction. The matching/orphan-deletion logic is unchanged from R1; only the table name and the new column write are different.

  Replace [pipeline.jl:473-533](../../packages/HimalayaUI/src/pipeline.jl) with:
  ```julia
  function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
      rows = Tables.rowtable(DBInterface.execute(db,
          """SELECT e.filename, x.id AS experiment_id
             FROM exposures e
             JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
             WHERE e.id = ?""", [exposure_id]))
      isempty(rows) && error("exposure $exposure_id not found")
      filename      = rows[1].filename
      experiment_id = rows[1].experiment_id

      cfg = config_from_db(db, experiment_id)
      pattern_filename = replace(cfg.integration_pattern, "{name}" => filename)
      dat_path = joinpath(analysis_dir, pattern_filename)
      isfile(dat_path) || error("dat file not found: $dat_path")

      q, I, σ = load_dat(dat_path)
      peaks_result = Himalaya.findpeaks(q, I, σ)

      diff_update_auto_peaks!(db, exposure_id, peaks_result, I)

      # Sharpness for `add` curations is sampled fresh inside effective_peaks
      # from the current trace — never stored on the curation row, never stale.
      eff = effective_peaks(db, exposure_id, q, I)
      candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
      group = auto_group(candidates)

      persist_analysis!(db, exposure_id, q, I, peaks_result, candidates, group, eff)
  end
  ```

- [ ] **Step 3: Rewrite `persist_analysis!` / `_persist_analysis_inner!`**

  The function signature gains `eff::NamedTuple` (the precomputed effective peaks). The body simplifies:
  - The excluded-q snapshot/restore block is gone (curations live in their own table, untouched).
  - Manual-peak ID merging is gone (curation-add ids are passed through `eff.peak_id` + `eff.peak_kind`).
  - The `q_to_peak_id` lookup in the candidate-persistence loop is replaced by zip-iterating `eff.q`, `eff.peak_id`, and `eff.peak_kind` together. When `Himalaya.indexpeaks` returns an `IndexEntry` whose `peaks` SparseVector contains a q-value found in `eff.q`, the corresponding `eff.peak_id[i]` and `eff.peak_kind[i]` go into the new `index_peaks` row's `peak_id` and `peak_kind` columns.

  `index_peaks.peak_id` references can be either an `auto_peaks.id` or a `peak_curations.id`; the `peak_kind` column (added in R2.1 step 2 — *not* re-added here) disambiguates. Future inserts populate `peak_kind` from the symbol carried in `eff.peak_kind`.

  This is the largest single edit in this task. Reference the current implementation closely; preserve every test that depends on speculative reattachment, custom-group reattachment, and stale-recovery behavior. Sketch of the surviving structure:

  ```julia
  function _persist_analysis_inner!(db, exposure_id, q, I,
                                     peaks_result, candidates, group, eff)
      # GONE: snapshot of `excluded_qs` from peaks(source='auto', excluded=1)
      # GONE: DELETE FROM peaks WHERE source='auto' + re-insert loop
      # GONE: re-application of was_excluded by q
      # GONE: manual peak ID merge (`peaks(source='manual')` no longer exists)
      # GONE: q_to_peak_id Dict built from peaks_result + manual rows

      # KEEPS: the speculative-reattachment branch — but its peak-id lookups
      # now consult `eff.peak_id` / `eff.peak_kind` (parallel arrays indexed
      # by position in `eff.q`) instead of querying `peaks` by source/q.
      # Build a small Dict for O(1) q→(peak_id, peak_kind) lookups inside
      # the candidate-persistence loop:
      eff_lookup = Dict{Float64, Tuple{Int, Symbol}}()
      for i in eachindex(eff.q)
          eff_lookup[eff.q[i]] = (eff.peak_id[i], eff.peak_kind[i])
      end

      # KEEPS: index_groups + index_group_members logic, custom-group
      # preservation, candidate persistence loop. The only change inside the
      # loop: when inserting an `index_peaks` row, look up via eff_lookup
      # and write both peak_id and peak_kind.

      # KEEPS: speculative re-resolution (rebuilding speculative `index_peaks`
      # rows whose anchor q still resolves to a peak in `eff`). The lookup
      # uses `eff_lookup` with tolerance matching SNAP_TOL.

      # KEEPS: status normalization (deleting orphan candidate rows that
      # don't survive auto_group / remove_subsets).
  end
  ```

  Roughly: ~180 lines of snapshot/restore/merge logic collapse into one Dict construction and a `peak_kind` write at one INSERT site. The remaining structure (groups, speculatives, status) is unchanged.

- [ ] **Step 4: Rewrite peak routes**

  - `GET /api/exposures/:id/peaks` ([routes_peaks.jl:19-26](../../packages/HimalayaUI/src/routes_peaks.jl)) — UNION ALL between `auto_peaks` and `peak_curations(kind='add')`, with synthesized `source` and `excluded` fields. Excluded auto peaks are still returned (with `excluded: true`) — the frontend's UI still shows them.

    ```julia
    @get "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
                   'auto' AS source,
                   CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
            FROM auto_peaks a
            LEFT JOIN peak_curations c
                ON c.exposure_id = a.exposure_id
               AND c.kind = 'exclude'
               AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
            WHERE a.exposure_id = ?
            UNION ALL
            SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence, NULL AS sharpness,
                   'manual' AS source, 0 AS excluded
            FROM peak_curations
            WHERE exposure_id = ? AND kind = 'add'
            ORDER BY q
        """, [id, id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows; bool_keys=(:excluded,))))
    end
    ```

  - `POST /api/exposures/:id/peaks` (manual add) — INSERT into `peak_curations` with `kind='add'`. `created_by` populated via `get_or_create_user!` from `X-Username`. Drop the `mark_all_indices_stale!` call (R3 will derive staleness from hash mismatch; until R3 lands, accept that other peoples' tabs may show non-stale indices for a moment after a peak edit — current `autoReanalyze` chain still fires).

  - `PATCH /api/peaks/:id { excluded: true }` — disambiguates by querying both tables for the id. If auto, INSERT into `peak_curations` with `kind='exclude'` (or DELETE if `excluded=false`). If curation-add, return 400 ("can't exclude a manually-added peak; delete it instead").

  - `DELETE /api/peaks/:id` — same disambiguation. If auto, refuse with 400 (auto peaks can only be excluded, not deleted). If curation-add, DELETE the `peak_curations` row.

- [ ] **Step 5: Update `json.jl` if needed for the new shape**

  Probably no change — `rows_to_json` already handles the joined shape via the existing `bool_keys` mechanism.

- [ ] **Step 6: Test the full curation lifecycle**

  In `test_routes_peaks.jl`, add testsets:
  - `POST /peaks` followed by `GET /peaks` returns the manual peak with `source='manual'`
  - `PATCH /peaks/:id { excluded: true }` on an auto peak inserts a curation row
  - `PATCH /peaks/:id { excluded: false }` on the same auto peak removes the curation row
  - Idempotence: two consecutive `PATCH { excluded: true }` produce one curation row, not two
  - `DELETE /peaks/:id` on a manual peak removes the curation row; subsequent `GET` doesn't return it
  - `DELETE /peaks/:id` on an auto peak returns 400

- [ ] **Step 7: Update existing pipeline testsets**

  The "incorporates manual peaks" and "ignores excluded auto peaks" testsets at [test_pipeline.jl:168, 210](../../packages/HimalayaUI/test/test_pipeline.jl) need their setup code updated to insert into `peak_curations` instead of `peaks(source='manual'|excluded=1)`. The behavioral assertions stay the same.

- [ ] **Step 8: Run full test suite**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  ```

- [ ] **Step 9: Commit**

  ```bash
  git add -p packages/HimalayaUI/src packages/HimalayaUI/test
  git commit -m "refactor(pipeline): cut analysis pipeline over to auto_peaks + peak_curations

  - analyze_exposure! reads effective_peaks (auto - excludes ∪ adds)
  - persist_analysis! threads effective peak ids through index_peaks
  - peak routes serve a joined view; frontend contract unchanged
  - manual-peak source-of-truth moves to peak_curations; auto-peak
    survival of curation events stays through diff-update from R1"
  ```

### Task R2.3: Drop legacy `peaks` table

This task only fires after R2.1 + R2.2 have been deployed and verified. It's separated for safety: ops can run R2.1+R2.2 in production, watch for issues, and run R2.3 once they're confident no rollback is needed.

**Migration impact:** The legacy `peaks` table is dropped permanently. Once this lands, restoring from a pre-R2 backup is the only path back.

- [ ] **Step 1: Verify R2.1+R2.2 have been live in the target environment for at least one sustained-use period**

  Operator-driven, not automated. The agent running this plan should pause here and ask the operator to confirm.

- [ ] **Step 2: Add a sentinel migration step**

  In `migrate_schema!`:
  ```julia
  # The DROP TABLE in migrate_r2_split_peaks! handles the cutover; this is a
  # belt-and-suspenders verifier that fails loudly if `peaks` is still around
  # after the migration was supposed to drop it.
  legacy = Tables.rowtable(DBInterface.execute(db,
      "SELECT 1 FROM sqlite_master WHERE type='table' AND name='peaks'"))
  if !isempty(legacy)
      @warn "Legacy 'peaks' table still present after R2 migration — investigate"
  end
  ```

  No code change beyond this verification — `migrate_r2_split_peaks!` already drops the table.

- [ ] **Step 3: Commit**

  ```bash
  git add packages/HimalayaUI/src/db.jl
  git commit -m "chore(db): post-R2 sentinel verifying legacy peaks table is dropped"
  ```

---

## R3: Content-hash memoization

**Migration impact:** Three columns added (`exposures.trace_hash`, `exposures.analysis_inputs_hash`, `indices.inputs_hash`). All start NULL on existing rows and populate on next analysis. The `mark_all_indices_stale!` mechanism is removed; staleness becomes a derived comparison.

### Task R3.1: SHA dependency + hash helpers

**Files:**
- Modify: `packages/HimalayaUI/Project.toml`
- Create: `packages/HimalayaUI/src/hash.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (include the new file)
- Create: `packages/HimalayaUI/test/test_hash.jl`

- [ ] **Step 1: Add SHA stdlib dep**

  In `Project.toml`:
  ```toml
  [deps]
  SHA = "ea8e919c-243c-51af-8825-aaa63cd721ce"
  ```

- [ ] **Step 2: Write the failing test**

  In `test_hash.jl`:
  ```julia
  using Test, HimalayaUI
  import HimalayaUI: hash_trace_file, hash_peak_set

  @testset "hash_trace_file is stable and content-addressed" begin
      mktempdir() do dir
          p = joinpath(dir, "t.dat")
          write(p, "0.1 1.0 0.1\n0.2 2.0 0.1\n")
          h1 = hash_trace_file(p)
          h2 = hash_trace_file(p)
          @test h1 == h2
          @test length(h1) == 64  # hex SHA-256

          # Different content → different hash
          write(p, "0.1 1.0 0.1\n0.2 99.0 0.1\n")
          h3 = hash_trace_file(p)
          @test h1 != h3
      end
  end

  @testset "hash_peak_set is order-independent and content-addressed" begin
      a = (q = [0.1, 0.2], sharpness = [1.0, 2.0])
      b = (q = [0.2, 0.1], sharpness = [2.0, 1.0])
      @test hash_peak_set(a) == hash_peak_set(b)

      c = (q = [0.1, 0.2], sharpness = [1.0, 2.0001])
      @test hash_peak_set(a) != hash_peak_set(c)
  end
  ```

- [ ] **Step 3: Implement**

  In `hash.jl`:
  ```julia
  using SHA

  """
      hash_trace_file(path) -> String

  SHA-256 over the bytes of the .dat file. Stable across Julia versions and
  re-reads of the same file. Used to detect "is the trace input unchanged
  since findpeaks last ran?"
  """
  function hash_trace_file(path::AbstractString)::String
      bytes2hex(open(SHA.sha256, path))
  end

  """
      hash_peak_set(eff::NamedTuple) -> String

  Content hash of an effective peak set. Inputs are sorted by q, encoded as
  consecutive (Float64, Float64) tuples in native byte order, and SHA-256ed.
  Used to detect "is the indexpeaks input unchanged since indices were
  computed?" Order-independent (sort key is q).
  """
  function hash_peak_set(eff::NamedTuple)::String
      n = length(eff.q)
      buf = Vector{UInt8}(undef, 16n)
      perm = sortperm(eff.q)
      for (i, k) in enumerate(perm)
          q  = Float64(eff.q[k])
          sh = Float64(eff.sharpness[k])
          unsafe_copyto!(buf, 16(i-1) + 1, reinterpret(UInt8, [q]), 1, 8)
          unsafe_copyto!(buf, 16(i-1) + 9, reinterpret(UInt8, [sh]), 1, 8)
      end
      bytes2hex(SHA.sha256(buf))
  end
  ```

- [ ] **Step 4: Wire into `HimalayaUI.jl`**

  Add `include("hash.jl")` in the module entry; export `hash_trace_file`, `hash_peak_set` (or leave unexported and access via `HimalayaUI.hash_*` in tests — match existing convention).

- [ ] **Step 5: Run tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  git add packages/HimalayaUI/Project.toml packages/HimalayaUI/Manifest.toml \
          packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/src/hash.jl \
          packages/HimalayaUI/test/test_hash.jl
  git commit -m "feat: hash helpers (SHA-256 trace + canonical peak-set hash)"
  ```

### Task R3.2: Schema columns + memoization in pipeline

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/src/pipeline.jl`
- Modify: `packages/HimalayaUI/src/routes_peaks.jl` (drop `mark_all_indices_stale!` calls)
- Modify: `packages/HimalayaUI/src/routes_analysis.jl` (return inputs_hash)
- Modify: `packages/HimalayaUI/test/test_pipeline.jl`

- [ ] **Step 1: Add columns + migration step**

  In `db.jl`'s `SCHEMA` constant, *edit the existing `CREATE TABLE exposures (...)` and `CREATE TABLE indices (...)` statements in place* — append the new columns to those table definitions. Do NOT add a separate ALTER-style block to `SCHEMA`; `SCHEMA` describes the post-migration state for fresh DBs (where the `migrate_schema!` ALTER calls below are no-ops because the columns are already present from CREATE).

  ```sql
  -- Inside the existing CREATE TABLE exposures (...):
  trace_hash           TEXT,
  analysis_inputs_hash TEXT,
  -- Inside the existing CREATE TABLE indices (...):
  inputs_hash          TEXT,
  ```

  In `migrate_schema!`:
  ```julia
  for stmt in [
      "ALTER TABLE exposures ADD COLUMN trace_hash TEXT",
      "ALTER TABLE exposures ADD COLUMN analysis_inputs_hash TEXT",
      "ALTER TABLE indices ADD COLUMN inputs_hash TEXT",
  ]
      try DBInterface.execute(db, stmt) catch end
  end
  ```

- [ ] **Step 2: Hash-guard `analyze_exposure!`**

  ```julia
  function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
      # ... resolve dat_path as before ...

      new_trace_hash = hash_trace_file(dat_path)
      stored_trace_hash = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT trace_hash FROM exposures WHERE id = ?", [exposure_id]))).trace_hash
      stored_trace_hash = ismissing(stored_trace_hash) ? nothing : String(stored_trace_hash)

      autopeaks_count = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT COUNT(*) AS n FROM auto_peaks WHERE exposure_id = ?", [exposure_id]))).n

      q, I, σ = load_dat(dat_path)
      if stored_trace_hash != new_trace_hash || autopeaks_count == 0
          peaks_result = Himalaya.findpeaks(q, I, σ)
          diff_update_auto_peaks!(db, exposure_id, peaks_result, I)
          DBInterface.execute(db,
              "UPDATE exposures SET trace_hash = ? WHERE id = ?",
              [new_trace_hash, exposure_id])
      end

      eff = effective_peaks(db, exposure_id, q, I)
      new_inputs_hash = hash_peak_set(eff)

      stored_inputs_hash = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT analysis_inputs_hash FROM exposures WHERE id = ?", [exposure_id]))).analysis_inputs_hash
      stored_inputs_hash = ismissing(stored_inputs_hash) ? nothing : String(stored_inputs_hash)

      indices_count = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT COUNT(*) AS n FROM indices WHERE exposure_id = ?", [exposure_id]))).n

      if stored_inputs_hash != new_inputs_hash || indices_count == 0
          # We need peaks_result here too for index_peaks linkage; if we skipped
          # findpeaks above, we still need to load it. The cheapest path: query
          # auto_peaks, synthesize a minimal NamedTuple equivalent.
          peaks_result_for_persist = synthesize_peaks_result(db, exposure_id, q, I)
          candidates = Himalaya.indexpeaks(eff.q, eff.sharpness)
          group = auto_group(candidates)
          persist_analysis!(db, exposure_id, q, I, peaks_result_for_persist,
                            candidates, group, eff)
          DBInterface.execute(db,
              "UPDATE exposures SET analysis_inputs_hash = ? WHERE id = ?",
              [new_inputs_hash, exposure_id])
          DBInterface.execute(db,
              "UPDATE indices SET inputs_hash = ? WHERE exposure_id = ?",
              [new_inputs_hash, exposure_id])
      end
  end
  ```

  `synthesize_peaks_result` reads `auto_peaks` rows for the exposure and reconstructs a NamedTuple matching `Himalaya.findpeaks`'s output shape:

  ```julia
  function synthesize_peaks_result(db::SQLite.DB, exposure_id::Int,
                                    q::Vector{Float64}, I::Vector{Float64})
      rows = Tables.rowtable(DBInterface.execute(db,
          """SELECT q, prominence, sharpness, findpeaks_index
             FROM auto_peaks WHERE exposure_id = ? ORDER BY q""", [exposure_id]))
      qs    = Float64[Float64(r.q)                                      for r in rows]
      proms = Float64[ismissing(r.prominence) ? 0.0 : Float64(r.prominence) for r in rows]
      shs   = Float64[ismissing(r.sharpness)  ? 0.0 : Float64(r.sharpness)  for r in rows]
      idxs  = Int[
          # Prefer the persisted findpeaks_index (exact local-maximum sample).
          # Fall back to nearest-grid-index for legacy rows from R2 migration
          # where findpeaks_index is NULL — these heal on next diff_update.
          ismissing(r.findpeaks_index) ?
              argmin(abs.(q .- Float64(r.q))) :
              Int(r.findpeaks_index)
          for r in rows
      ]
      (q = qs, indices = idxs, prominence = proms, sharpness = shs)
  end
  ```

  The persisted `findpeaks_index` keeps the reconstruction bit-equivalent to a fresh `findpeaks` call (the local-maximum grid sample, not its argmin-nearest neighbour, which can differ for asymmetric peaks). `diff_update_auto_peaks!` (R1) is the writer — it must set `findpeaks_index` whenever it inserts or updates an auto-peak row.

  Field shape comes from [src/peakfinding.jl:100-103](../../src/peakfinding.jl) — `findpeaks` returns `(indices, q, prominence, sharpness)` as a NamedTuple. The `synthesize_peaks_result` reconstruction matches that shape exactly.

- [ ] **Step 3: Drop `mark_all_indices_stale!` and clean up orphan `'stale'` status values**

  Delete the function from [routes_peaks.jl:9-16](../../packages/HimalayaUI/src/routes_peaks.jl) and all three call sites (lines 38, 78, 117). Frontend banner gates on hash mismatch in the next step.

  R2's migration set `status='stale'` on every existing index to force reanalysis under the new effective_peaks model. Once R3 derives staleness from hash mismatch, those `'stale'` values are dead — the column still holds `'candidate'` for live indices, but rows from the R2 transition would otherwise carry an orphaned value forever. As part of R3.2's `migrate_r3_hash_columns!` (the new column-add migration step), normalize:
  ```julia
  # Once analysis_inputs_hash exists and indices.inputs_hash is the source of
  # truth for staleness, the old 'stale' enum value is dead. Any existing
  # 'stale' rows came from R2's transitional UPDATE; normalize them to
  # 'candidate' (the next analyze run will set inputs_hash; until then,
  # hash mismatch with NULL on the index renders them as stale to the UI).
  DBInterface.execute(db, "UPDATE indices SET status = 'candidate' WHERE status = 'stale'")
  ```

  The `status` column itself stays — it's used for `'candidate'` and may grow other values later. We're just retiring the `'stale'` value.

- [ ] **Step 4: Surface `analysis_inputs_hash` on responses**

  - `GET /api/exposures/:id`: include `analysis_inputs_hash` in the JSON.
  - `GET /api/exposures/:id/indices`: include the exposure's current `analysis_inputs_hash` and each index's `inputs_hash` so the frontend can compare.
  - `GET /api/samples/:id/exposures`: include `analysis_inputs_hash` per row.

- [ ] **Step 5: Frontend banner update**

  In `StaleIndicesBanner.tsx`, change the condition from `status === "stale"` to:
  ```typescript
  const isStale = exposure?.analysis_inputs_hash !== undefined
                && indices.some(ix => ix.inputs_hash !== exposure.analysis_inputs_hash);
  ```

  Update Vitest tests accordingly. The banner UI behavior (Re-analyze button posting to `/api/exposures/:id/analyze`) stays unchanged.

- [ ] **Step 6: Add memoization regression tests**

  In `test_pipeline.jl`:
  ```julia
  @testset "analyze_exposure! skips findpeaks when trace_hash unchanged" begin
      # Run analyze once, instrument findpeaks call count, run again, assert
      # the count didn't increment. Easiest: monkey-patch via Mocking.jl, or
      # add a counter to a test-only module-level Ref that findpeaks bumps.
      # If neither is convenient, assert by elapsed time (re-run is dramatically
      # faster than first run on a real fixture).
  end

  @testset "analyze_exposure! skips indexpeaks when peak set hash unchanged" begin
      # Similar pattern.
  end

  @testset "analyze_exposure! re-runs findpeaks when trace bytes change" begin
      # Touch the .dat file in-between, assert findpeaks ran again.
  end
  ```

- [ ] **Step 7: Run tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  cd packages/HimalayaUI/frontend && npm test && npm run build && cd -
  git add -p
  git commit -m "feat: hash-guarded analysis memoization

  - exposures.trace_hash skips findpeaks when .dat unchanged
  - exposures.analysis_inputs_hash skips indexpeaks when effective peak set unchanged
  - indices.inputs_hash drives StaleIndicesBanner via hash mismatch
  - mark_all_indices_stale! removed (staleness is now derived)"
  ```

---

## R4: Promote `user_actions` to `curation_events`

**Migration impact:** `user_actions` gains two columns (`payload`, `undoes_event_id`) and one index. Historical rows have `payload = NULL` — those events can't replay deterministically, but they're never used as the source of truth (the materialized views are pre-populated by R2 and remain authoritative for state up to this point). New events from this PR forward have full payloads.

The materialized-view contract is the new load-bearing rule: every curation route MUST go through `apply_event!` and `apply_event!` MUST update both the log and the views in one transaction.

### Task R4.1: Schema + `apply_event!` helper

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Create: `packages/HimalayaUI/src/events.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/src/actions.jl`
- Create: `packages/HimalayaUI/test/test_events.jl`

- [ ] **Step 1: Schema migration**

  ```julia
  # in migrate_schema!
  for stmt in [
      "ALTER TABLE user_actions ADD COLUMN payload TEXT",
      "ALTER TABLE user_actions ADD COLUMN undoes_event_id INTEGER REFERENCES user_actions(id)",
  ]
      try DBInterface.execute(db, stmt) catch end
  end
  ```

  Add to `SCHEMA`:
  ```sql
  CREATE INDEX IF NOT EXISTS idx_events_by_exposure
      ON user_actions(entity_type, entity_id, id);
  ```

- [ ] **Step 2: `apply_event!` helper**

  In `events.jl`:
  ```julia
  using JSON3, SQLite, DBInterface

  """
      apply_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id=nothing)
        -> event_id::Int

  Atomic event-append + view-update. The log and the views must move together
  or neither moves. Returns the newly-inserted event id.

  `payload` is any JSON-serializable Dict / NamedTuple / nothing. If nothing,
  the event is recorded but no view update fires (use sparingly — most actions
  should carry a payload).
  """
  function apply_event!(db::SQLite.DB, req;
                        kind::String,
                        entity_type::String,
                        entity_id::Integer,
                        payload = nothing,
                        undoes_event_id::Union{Int,Nothing} = nothing)
      username = get_username(req)
      user_id  = username === nothing ? nothing : get_or_create_user!(db, username)
      payload_json = payload === nothing ? nothing : JSON3.write(payload)

      event_id = SQLite.transaction(db) do
          res = DBInterface.execute(db,
              """INSERT INTO user_actions
                 (user_id, action, entity_type, entity_id, payload, undoes_event_id)
                 VALUES (?, ?, ?, ?, ?, ?)""",
              [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id])
          eid = Int(DBInterface.lastrowid(res))
          payload === nothing || update_view_for_event!(db, kind, entity_id, payload, eid)
          eid
      end

      # Best-effort SSE broadcast — fires AFTER the transaction commits, so a
      # subscriber never sees an event that was rolled back. If the process
      # dies between commit and broadcast, the event is durable in user_actions
      # but the frame is lost; clients reconcile on reconnect via TanStack
      # Query refetch (see R5a). isdefined check keeps the helper a soft
      # dependency: R4 lands before R5a, and apply_event! must work without
      # broadcast wired up yet.
      if isdefined(@__MODULE__, :broadcast_event!)
          try
              broadcast_event!(event_id, kind, entity_type, Int(entity_id), user_id, payload_json)
          catch err
              @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
          end
      end
      event_id
  end

  """
  Dispatcher that updates materialized views based on event kind.
  Currently a no-op for most events; populated as routes migrate to apply_event!.

  **Payload contract:** `payload` is normalized to a `JSON3.Object` (or `nothing`)
  before this is called. The live path in `apply_event!` round-trips through
  `JSON3.write` → `JSON3.read` so dispatcher branches see the same shape they'd
  see during `rebuild_views_from_log!`. This eliminates the Symbol-key vs
  String-key footgun: `JSON3.Object` supports both `obj.q` and `obj[:q]` /
  `obj["q"]`. Branch code accesses fields uniformly via `payload.q` style.
  """
  function update_view_for_event!(db, kind, entity_id, payload, event_id)
      # Implementations added per route in subsequent tasks. Initially:
      kind == "noop_test" && return
      # default: no view update
  end
  ```

  Update `apply_event!` to canonicalize before dispatch:
  ```julia
  # Inside the SQLite.transaction block, before the dispatcher call:
  payload_canonical = payload_json === nothing ? nothing : JSON3.read(payload_json)
  payload_canonical === nothing || update_view_for_event!(db, kind, entity_id, payload_canonical, eid)
  ```

  Replace the earlier line `payload === nothing || update_view_for_event!(db, kind, entity_id, payload, eid)` with the canonicalized version. This costs one parse per event on the live path; the upside is `rebuild_views_from_log!` and the live dispatcher follow exactly the same code path, so the property test isn't fooled by type-shape differences.

- [ ] **Step 3: `log_action!` becomes a thin wrapper**

  In `actions.jl`, leave `log_action!` in place but redirect:
  ```julia
  function log_action!(db, req; action, entity_type, entity_id, note=nothing)
      apply_event!(db, req;
                   kind = action,
                   entity_type = entity_type,
                   entity_id = entity_id,
                   payload = note === nothing ? nothing : Dict(:note => note))
      nothing
  end
  ```

  Existing call sites continue to work; new code can call `apply_event!` directly with structured payloads.

- [ ] **Step 4: Test the round-trip**

  In `test_events.jl`:
  ```julia
  @testset "apply_event! writes log and runs view update atomically" begin
      mktempdir() do dir
          db = HimalayaUI.open_db(joinpath(dir, "h.db"))
          # HTTP.Request constructor is (method, target, headers, body); pass
          # headers as Vector{Pair} so HTTP.header(req, "X-Username") resolves.
          req = HTTP.Request("POST", "/x",
              ["X-Username" => "alice"], UInt8[])

          eid = HimalayaUI.apply_event!(db, req;
              kind = "test_kind", entity_type = "exposure", entity_id = 42,
              payload = Dict(:foo => "bar"))
          @test eid > 0

          row = first(Tables.rowtable(DBInterface.execute(db,
              "SELECT * FROM user_actions WHERE id = ?", [eid])))
          @test String(row.action) == "test_kind"
          @test JSON3.read(row.payload).foo == "bar"
      end
  end

  @testset "log_action! still works (legacy wrapper)" begin
      # Backwards compat: existing call sites unchanged.
  end
  ```

- [ ] **Step 5: Run tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  git add -p
  git commit -m "feat(events): apply_event! infrastructure with atomic view updates

  Promotes user_actions to a structured event log. log_action! becomes a thin
  wrapper for backwards compat; new code carries structured payloads."
  ```

### Task R4.2: Migrate per-route call sites with payloads

For each of the ~20 `log_action!` call sites, decide what payload makes the resulting event useful for replay and instrumentation. Reference the spec's payload-discipline table — recording a field at design time is cheap; reconstructing it later is impossible.

**Migration impact:** No schema change. Routes change shape: they validate input and call `apply_event!` *with no view-side INSERT/UPDATE/DELETE*. The view-update dispatcher in `events.jl` is the sole writer to view tables. This is the contract enforcer for `rebuild_views_from_log!`: re-folding the log produces the same view state because the writes only ever happened through the dispatcher in the first place.

**Atomicity requirement (per kind):** for each event kind, the route refactor and the matching dispatcher branch MUST land in the same commit. The danger is the in-between state: a route flipped to `apply_event!` *before* its dispatcher branch lands logs the event, returns 200 to the client, and silently drops the view-side write — pure data loss inside R4.2. Do not split "convert routes" from "implement dispatcher branches" into two passes. The `rebuild_views_from_log!` property test only catches this if it covers every kind, which is exactly why it's the per-kind regression.

This task is mechanical but high-volume. Each route gets its own commit (route + dispatcher branch + property test together).

For each route in:
- `routes_peaks.jl`: `add_peak`, `exclude_peak`/`include_peak`, `remove_peak`, `delete_auto_peak_legacy`
- `routes_analysis.jl`: `confirm_index`, `exclude_index`, `create_speculative`, `delete_speculative`
- `routes_exposures.jl`: `set_status`, `select_exposure`, `add_tag`, `remove_tag`, `analyze`
- `routes_samples.jl`: `update_sample`, `add_tag`, `remove_tag`
- `routes_messages.jl`: `add_message`
- `routes_experiments.jl`: `update_experiment`, `analyze`, `reingest`

- [ ] **Step (per route): refactor route + add dispatcher branch**

  Routes get strictly *smaller* — they validate input, then call `apply_event!`. They do NOT INSERT into `peak_curations` or any other view-mat table directly.

  For example, `POST /api/exposures/:id/peaks` (add manual peak):
  ```julia
  @post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
      db   = current_db()
      body = json(req)
      q    = Float64(body.q)

      # Validate inputs (could check exposure exists, q is finite, etc.).
      # Then delegate everything to apply_event!.
      event_id = apply_event!(db, req;
          kind = "peak_added",
          entity_type = "exposure",
          entity_id = id,
          payload = Dict(:q => q))

      # Return the new curation row by querying after dispatch.
      cur = first(Tables.rowtable(DBInterface.execute(db,
          """SELECT id, exposure_id, q FROM peak_curations
             WHERE exposure_id = ? AND kind = 'add'
             ORDER BY id DESC LIMIT 1""", [id])))
      HTTP.Response(201, ["Content-Type" => "application/json"],
          JSON3.write(Dict(:id => Int(cur.id), :exposure_id => id,
                            :q => Float64(cur.q), :source => "manual",
                            :event_id => event_id)))
  end
  ```

  And the corresponding dispatcher branch in `update_view_for_event!`:
  ```julia
  if kind == "peak_added"
      # payload is a JSON3.Object (see canonicalization in apply_event!).
      DBInterface.execute(db,
          """INSERT INTO peak_curations (exposure_id, kind, q, created_by)
             VALUES (?, 'add', ?, (SELECT user_id FROM user_actions WHERE id = ?))""",
          [Int(entity_id), Float64(payload.q), event_id])
      return
  end

  # Second worked example: peak_unexcluded shows the undoes_event_id pattern.
  # The route reads the original peak_excluded event id, passes it as
  # undoes_event_id to apply_event!, and the dispatcher uses that plus the
  # payload to remove the matching exclude curation.
  if kind == "peak_unexcluded"
      # payload.q is snapshotted from the original peak_excluded event (looked
      # up via undoes_event_id by the route), so on the same trace the value
      # round-trips bit-equal through JSON3 and exact equality is sufficient.
      # Guard with a tiny tolerance anyway — same shape as effective_peaks's
      # is_excluded — so a manually-crafted replay/test can't fail for
      # rounding reasons. Cheap insurance, not load-bearing.
      DBInterface.execute(db,
          """DELETE FROM peak_curations
             WHERE exposure_id = ? AND kind = 'exclude'
               AND ABS(q - ?) <= MAX(1e-6, ABS(q) * 0.001)""",
          [Int(entity_id), Float64(payload.q)])
      return
  end
  ```

  Note the entity_type change from `"peak"` to `"exposure"` — folding the log per exposure is the common access pattern. Don't change `entity_type` for events that legitimately don't belong to an exposure (e.g., `update_experiment`).

  The `(SELECT user_id FROM user_actions WHERE id = ?)` self-reference depends on `event_id` being the just-inserted row's id — guaranteed because `apply_event!` calls the dispatcher *inside* the same transaction, after the INSERT. If the dispatcher is ever invoked outside `apply_event!` (e.g., by `rebuild_views_from_log!`), the row already exists in `user_actions` so the lookup still resolves. Adding an explicit guard is unnecessary; calling the dispatcher with a stale or fabricated `event_id` is a programmer error worth letting fail loudly via NULL `created_by` (which becomes a visible "unattributed" row in the UI).

  **Migration of `entity_type` for historical rows:** the `idx_events_by_exposure` index on `(entity_type, entity_id, id)` is only useful if peak/index events all key under `entity_type='exposure'`. Pre-R4 rows logged via `log_action!` use `entity_type='peak'` (in `routes_peaks.jl` for add/exclude/include/remove) and `entity_type='index'` (in `routes_analysis.jl` for confirm/exclude/speculative). These need a one-shot rewrite during R4.1's migration step:

  ```julia
  # In migrate_schema! (R4.1), AFTER the payload/undoes_event_id ALTERs and
  # BEFORE creating idx_events_by_exposure. Rewrites historical rows so the
  # new index is useful for fold-by-exposure queries.
  function migrate_r4_rebase_entity_type!(db::SQLite.DB)
      # Sentinel: skip if we've already rebased (presence of peak/index rows
      # under their old entity_type is the signal).
      legacy = first(Tables.rowtable(DBInterface.execute(db, """
          SELECT COUNT(*) AS n FROM user_actions
          WHERE entity_type IN ('peak', 'index')
      """))).n
      Int(legacy) == 0 && return

      SQLite.transaction(db) do
          # Resolve peak → exposure via the legacy peaks table if it still
          # exists, else via the new auto_peaks/peak_curations split.
          peaks_exists = !isempty(Tables.rowtable(DBInterface.execute(db,
              "SELECT 1 FROM sqlite_master WHERE type='table' AND name='peaks'")))
          if peaks_exists
              DBInterface.execute(db, """
                  UPDATE user_actions SET
                      entity_type = 'exposure',
                      entity_id = (SELECT exposure_id FROM peaks WHERE peaks.id = user_actions.entity_id)
                  WHERE entity_type = 'peak'
                    AND EXISTS (SELECT 1 FROM peaks WHERE peaks.id = user_actions.entity_id)
              """)
          else
              DBInterface.execute(db, """
                  UPDATE user_actions SET
                      entity_type = 'exposure',
                      entity_id = COALESCE(
                          (SELECT exposure_id FROM auto_peaks      WHERE auto_peaks.id      = user_actions.entity_id),
                          (SELECT exposure_id FROM peak_curations  WHERE peak_curations.id  = user_actions.entity_id)
                      )
                  WHERE entity_type = 'peak'
                    AND EXISTS (
                        SELECT 1 FROM auto_peaks     WHERE auto_peaks.id     = user_actions.entity_id
                        UNION SELECT 1 FROM peak_curations WHERE peak_curations.id = user_actions.entity_id
                    )
              """)
          end
          # Index events: resolve via indices.exposure_id.
          DBInterface.execute(db, """
              UPDATE user_actions SET
                  entity_type = 'exposure',
                  entity_id = (SELECT exposure_id FROM indices WHERE indices.id = user_actions.entity_id)
              WHERE entity_type = 'index'
                AND EXISTS (SELECT 1 FROM indices WHERE indices.id = user_actions.entity_id)
          """)
          # Stragglers (peak/index id no longer resolves) keep their original
          # entity_type — they won't appear in fold-by-exposure queries but
          # they're not lost; they're queryable by raw kind/id.
      end
  end
  ```

  Wire this into `migrate_schema!` after the payload ALTERs and before the index CREATE so rebased rows land under the index correctly.

  The `created_by` is resolved by reading back from the just-inserted `user_actions` row — keeps the dispatcher self-contained without threading additional args through.

- [ ] **Step (test): assert `rebuild_views_from_log!` matches incremental state**

  After each route migrates, add a property test:
  ```julia
  @testset "rebuild_views_from_log! reproduces incremental state for $kind" begin
      # 1. Set up a fresh DB.
      # 2. Apply a sequence of events via the route handlers.
      # 3. Snapshot the view state (peak_curations rows, index_group_members rows).
      # 4. Drop the view tables, re-create empty.
      # 5. Run rebuild_views_from_log!(db, exposure_id).
      # 6. Assert view state is bit-identical to the snapshot.
  end
  ```

  This test is the contract enforcer. **If a route still writes a view directly (bypasses `apply_event!`), this test fails** because the rebuild won't reproduce the row that the route inserted.

- [ ] **Step: speculative index events**

  `create_speculative`'s payload should record the full ratio_to_peak_id mapping, the resulting basis/score/r²/lattice_d, AND the IDs and scores of competing auto candidates at the moment of creation. This supports the "did users prefer the highest-scoring candidate?" instrumentation question.

- [ ] **Step: index_confirmed events**

  Critical payload field: array of competing candidate scores at confirmation time:
  ```julia
  payload = Dict(
      :index_id => index_id,
      :phase => phase_name,
      :basis => basis,
      :score => score,
      :r_squared => r_squared,
      :competing => [
          Dict(:phase => p.phase, :basis => p.basis, :score => p.score, :r_squared => p.r_squared)
          for p in competing_candidates
      ],
  )
  ```

  Without `competing`, you lose the ability to ask "did users override the engine's top choice?" later.

- [ ] **Step: analyze_run event from pipeline.jl**

  Emit an `analyze_run` event **unconditionally** at the end of `analyze_exposure!` (whether or not findpeaks/indexpeaks ran). The point is observability of the memoization layer — a no-op run that hits both hash caches is exactly the data point that proves R3 is doing its job.

  Wrap the body to capture timing:
  ```julia
  function analyze_exposure!(db::SQLite.DB, exposure_id::Int, analysis_dir::String)
      t0 = time()
      # ... existing body: hash check, diff_update, effective_peaks, persist ...

      # The Boolean fields use the values captured BEFORE the conditional
      # blocks ran, so a hit/miss is recorded faithfully.
      duration_ms = round(Int, (time() - t0) * 1000)
      apply_event!(db, _system_request();
          kind = "analyze_run",
          entity_type = "exposure",
          entity_id = exposure_id,
          payload = Dict(
              :trace_hash_before  => stored_trace_hash,
              :trace_hash_after   => new_trace_hash,
              :inputs_hash_before => stored_inputs_hash,
              :inputs_hash_after  => new_inputs_hash,
              :findpeaks_skipped  => (stored_trace_hash  == new_trace_hash),
              :indexpeaks_skipped => (stored_inputs_hash == new_inputs_hash),
              :duration_ms        => duration_ms,
              :effective_peaks_count => length(eff.q),  # not "auto" — eff
                                                         # excludes excluded
                                                         # autos and includes
                                                         # add curations.
          ))
  end
  ```

  `_system_request()` returns a synthetic `HTTP.Request` with no `X-Username` so the event's `user_id` is NULL — pipeline runs aren't attributed to any actor:
  ```julia
  # In events.jl alongside apply_event!.
  _system_request() = HTTP.Request("INTERNAL", "/_system", Pair{String,String}[], UInt8[])
  ```
  `get_username(req)` returns `nothing` for a request with no `X-Username` header (see [actions.jl](../../packages/HimalayaUI/src/actions.jl)), which is the existing convention. (If a route-triggered analyze should be attributed, the route can call `apply_event!` itself with its own `req` for an `analyze_requested` event before invoking `analyze_exposure!`.)

  Rename `:auto_peaks_count` → `:effective_peaks_count` since the value is `length(eff.q)`, which counts post-exclusion + add curations. The earlier name was misleading.

  This payload supports both pipeline instrumentation (memoization hit rate, analysis duration distribution) and quality metrics (effective-peak count over time).

- [ ] **Step: rebuild-from-log helper**

  Add `rebuild_views_from_log!(db, exposure_id)` to `events.jl`. Folds the log to reconstruct `peak_curations` and `index_group_members` from scratch. Used for migration and disaster recovery; tested against the property "after `apply_event!` of any sequence, the views match what `rebuild_views_from_log!` would produce."

  This testset becomes load-bearing: every event kind that updates a view must round-trip through rebuild correctly. If a payload is missing fields needed for replay, this test catches it.

- [ ] **Step: commit per route group**

  Group commits by route file:
  ```bash
  git add -p packages/HimalayaUI/src/routes_peaks.jl packages/HimalayaUI/src/events.jl
  git commit -m "feat(events): structured payloads for peak curation events"
  # ... one commit per routes_*.jl ...
  ```

---

## R5a: SSE broadcast + LWW everywhere

**Migration impact:** No schema. New endpoint, new client subscriber. Best-effort delivery (event is durable in `user_actions`; SSE frame may be lost if process dies between commit and broadcast). Clients reconcile on reconnect via TanStack Query refetch.

**Scope split:** R5a (this section) ships SSE broadcast and LWW conflict resolution everywhere. R5b adds `If-Match` + 409 retry on delta-shaped routes — gated on R4 instrumentation showing actual conflicts in practice. R5b may never ship if R4 data shows contention is rare.

### Task R5.0: Oxygen.jl streaming + heartbeat-timer spike

Before any R5a code commits, the spike must resolve **two** open questions:

1. **How to express SSE in Oxygen.jl 1.10.x.** The `@stream` macro mentioned in earlier drafts may not exist in this version.
2. **How to interleave event delivery with heartbeat frames.** `Threads.wait(::Condition)` doesn't natively support a timeout; a naive subscriber loop either waits forever (no heartbeats) or polls (defeats the point of a Condition). The right pattern (`Timer` that periodically calls `notify(SSE_WAKEUP)`, or a `take!(channel)` against a `Channel` that the heartbeat loop also pushes to) needs a working POC before R5.1 commits.

Both questions are coupled — the chosen Oxygen API affects what abstractions are available for the loop. Resolve together.

**Files:** Throwaway scratch script under `scratch/sse_spike.jl` (gitignored — see CLAUDE.md "Code layout" → `scratch/`).

- [ ] **Step 1: Read Oxygen.jl 1.10 docs / source for streaming**

  Determine the actual API for streaming responses. Likely candidates:
  - A `streamhandler` decorator
  - `@get` with manual `HTTP.Stream` access via `req.stream` or similar
  - Falling back to dropping into `HTTP.serve` for the SSE endpoint while keeping Oxygen for everything else

- [ ] **Step 2: Write a spike that exercises both delivery paths**

  Goal: a `/api/sse_test` endpoint that:
  - Holds the connection open against `curl -N`.
  - Receives notifications from a separate task (mimicking `broadcast_event!`) and emits an SSE event frame within 100 ms of the notify call.
  - Emits a `:heartbeat\n\n` comment frame every 15 s of silence.
  - Cleanly removes the subscriber on connection close (verify via `:strace`-style log of `eof(stream)` returning true).

  Verify both paths in one run: `curl -N` for ~40 s with no events (heartbeats fire), then a Julia REPL pushes 3 notifies in quick succession (events fire within 100 ms each), then close the curl (subscriber list shrinks).

- [ ] **Step 3: Decision**

  Document both decisions in `docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md` (one paragraph in §R5a Server side):
  1. The chosen Oxygen API (or fallback to separate `HTTP.serve`).
  2. The chosen wakeup mechanism (Timer-driven notify, Channel-based merge, etc.).

  If either question can't be resolved cleanly, do not proceed to R5.1 — instead surface the blocker before architecting code around assumptions that won't hold.

### Task R5.1: Server-side SSE infrastructure

**Files:**
- Modify: `packages/HimalayaUI/src/server.jl`
- Modify: `packages/HimalayaUI/src/events.jl` (broadcast on commit)

- [ ] **Step 1: Subscriber registry**

  In `server.jl`, add module-level state:
  ```julia
  # Each entry is (stream, pending::Channel{String}). The Channel is the
  # subscriber's single source of frames — both real events (pushed by
  # broadcast_event!) and heartbeats (pushed by a per-subscriber Timer set up
  # in Step 2's handler). The handler loop is just `for frame in pending`,
  # which blocks correctly without TOCTOU.
  const SSE_SUBSCRIBERS = Ref{Vector{Any}}([])
  const SSE_LOCK        = ReentrantLock()
  ```

  No shared `Condition` is needed: each subscriber blocks on its own `pending::Channel`, and broadcasting `put!`s onto every subscriber's channel directly. Earlier drafts used a shared `Threads.Condition` with `isready` + `wait`, which races (a notify between `isready` returning false and `wait` registering can be missed, blocking up to a heartbeat interval). Event latency is now bounded by network RTT, with no busy-poll and no race.

  *Conditional fallback:* if the R5.0 spike finds the Channel pattern doesn't compose with the chosen Oxygen streaming primitive, fall back to a shared `SSE_WAKEUP::Threads.Condition` plus per-subscriber `pending::Channel`, with the loop using `take!` inside the handler and `notify` after `put!` in `broadcast_event!`. The Channel-only design is the default; the Condition design is the documented escape hatch.

- [ ] **Step 2: SSE endpoint (using whichever API the spike landed on)**

  Pseudocode using a generic stream type — replace with the concrete API from Task R5.0:
  ```julia
  function sse_handler(stream)
      # Reverse-proxy-friendly headers. X-Accel-Buffering: no defeats nginx
      # buffering (otherwise frames stall until the response closes). The
      # heartbeat below keeps the connection alive past intermediaries' idle
      # timeouts (typically 60s).
      HTTP.setheader(stream, "Content-Type"      => "text/event-stream")
      HTTP.setheader(stream, "Cache-Control"     => "no-cache")
      HTTP.setheader(stream, "Connection"        => "keep-alive")
      HTTP.setheader(stream, "X-Accel-Buffering" => "no")
      HTTP.startwrite(stream)

      pending = Channel{String}(64)  # per-subscriber buffer
      lock(SSE_LOCK) do
          push!(SSE_SUBSCRIBERS[], (stream=stream, pending=pending))
      end

      # Heartbeat producer: an independent Timer that put!s a comment frame onto
      # `pending` every 15s. This collapses event delivery and heartbeats into a
      # single source — `take!(pending)` blocks correctly with no TOCTOU between
      # an `isready` check and `wait`. (Earlier drafts used `isready` + a shared
      # Condition, which races: a notify between `isready` returning false and
      # `wait` registering can be missed, blocking up to a heartbeat interval.)
      heartbeat_timer = Timer(15; interval = 15) do _
          try put!(pending, ":heartbeat\n\n") catch end  # closed channel: noop
      end

      try
          while !eof(stream)
              frame = take!(pending)  # blocks until event or heartbeat
              write(stream, frame)
          end
      finally
          close(heartbeat_timer)
          lock(SSE_LOCK) do
              filter!(sub -> sub.stream !== stream, SSE_SUBSCRIBERS[])
          end
          close(pending)
      end
  end
  ```

  Subscriber state matches Step 1: each subscriber has its own `pending::Channel`; `broadcast_event!` `put!`s onto every channel directly; no shared Condition is involved. Task R5.0 should validate this pattern composes with the chosen Oxygen streaming primitive; if not, fall back to the Condition variant documented in Step 1.

- [ ] **Step 3: Broadcast helper**

  In `events.jl`, after the transaction in `apply_event!` commits:
  ```julia
  broadcast_event!(event_id, kind, entity_type, entity_id, user_id, payload_json)
  ```

  Implement `broadcast_event!` to enqueue a frame on each subscriber's `pending` channel and notify the condition. Closed/full channels signal a dead subscriber that gets pruned on next iteration:
  ```julia
  # Trivial helper — kept here so the broadcast snippet is self-contained.
  function lookup_username(db::SQLite.DB, user_id::Integer)::Union{String, Nothing}
      rows = Tables.rowtable(DBInterface.execute(db,
          "SELECT name FROM users WHERE id = ?", [Int(user_id)]))
      isempty(rows) || ismissing(rows[1].name) ? nothing : String(rows[1].name)
  end

  # Frontend self-echo filter (R5a) compares `event.actor === username`. For
  # anonymous edits (no X-Username), both sides resolve to nullish, so every
  # anonymous user filters every other anonymous user's events. This is
  # acceptable under the X-Username trust model — anonymous users aren't
  # identified anyway. If lab usage ever introduces multi-tab anonymous
  # collaboration, switch to per-tab client ids; not worth the cost today.

  function broadcast_event!(event_id, kind, entity_type, entity_id, user_id, payload_json)
      actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
      msg = JSON3.write(Dict(
          :id => event_id, :kind => kind,
          :entity_type => entity_type, :entity_id => entity_id, :actor => actor,
          :payload => payload_json === nothing ? nothing : JSON3.read(payload_json),
      ))
      frame = "event: curation\ndata: $msg\n\n"
      lock(SSE_LOCK) do
          to_drop = []
          for sub in SSE_SUBSCRIBERS[]
              try
                  put!(sub.pending, frame)
              catch
                  push!(to_drop, sub)
              end
          end
          for sub in to_drop
              filter!(x -> x !== sub, SSE_SUBSCRIBERS[])
          end
      end
      # No notify: each subscriber blocks on its own pending Channel, and put!
      # above wakes its take! directly.
  end
  ```

  **Best-effort:** If the process dies between `apply_event!`'s `SQLite.transaction` commit and `broadcast_event!`, the event is durable but no frame is sent. Clients reconcile on reconnect via TanStack Query refetch. The principled fix (`Last-Event-ID` header → server-side replay from `user_actions`) is deferred to a v2 — see Open Questions in the spec.

- [ ] **Step 4: Tests**

  Add `test_sse.jl` covering:
  - Subscriber receives broadcast on `apply_event!` within 100ms (not 1s+).
  - Multiple subscribers each receive.
  - Closed streams are pruned from registry.
  - Heartbeat frame fires after 15s of silence.
  - `X-Accel-Buffering: no` header is set on the response.

### Task R5.2 (deferred to R5b): `If-Match` conflict checks

This task only fires after R4 instrumentation has been live for a sustained-use period (≥4 weeks of typical lab use, or ≥500 curation events across ≥10 distinct exposures, whichever lands first) and the analysis below crosses the threshold:

**Concrete go/no-go threshold:** Run the canonical query below against `user_actions`. It enumerates delta-shaped events (`peak_added`, `peak_excluded`, `peak_unexcluded`, `index_confirmed`, `index_unconfirmed`, `speculative_created`) and counts those with *at least one* prior delta-shaped event on the same exposure within 5 seconds, from a *different* actor. If `contended / total >= 0.02`, ship R5b. Otherwise, LWW under R5a is sufficient and R5b can be skipped indefinitely.

```sql
WITH delta_events AS (
    SELECT id, user_id, entity_id AS exposure_id, created_at
    FROM user_actions
    WHERE action IN (
        'peak_added', 'peak_excluded', 'peak_unexcluded',
        'index_confirmed', 'index_unconfirmed', 'speculative_created'
    )
      AND entity_type = 'exposure'  -- post-R4 entity_type rebase
      AND user_id IS NOT NULL       -- anonymous edits don't count as contention
),
contended AS (
    SELECT e.id
    FROM delta_events e
    WHERE EXISTS (
        SELECT 1 FROM delta_events p
        WHERE p.exposure_id = e.exposure_id
          AND p.user_id != e.user_id
          AND p.created_at >= datetime(e.created_at, '-5 seconds')
          AND p.created_at <  e.created_at
    )
)
SELECT
    (SELECT COUNT(*) FROM delta_events) AS total,
    (SELECT COUNT(*) FROM contended)    AS contended_count,
    CAST((SELECT COUNT(*) FROM contended) AS REAL)
        / NULLIF((SELECT COUNT(*) FROM delta_events), 0) AS rate;
```

Decision rule: `rate >= 0.02 AND total >= 500 AND span >= 4 weeks` → ship R5b. The semantic choices are pre-committed here so the threshold can't be reverse-engineered:
- "At least one" prior cross-actor event in the window — not "exactly one." Three-way concurrent edits count once per affected event, which is the right behavior (each is contended).
- "Different actor" excludes self-edits — a user editing their own work in rapid succession isn't multiplayer contention.
- Anonymous edits (`user_id IS NULL`) are excluded from both numerator and denominator — they can't be attributed and would inflate the denominator without ever being counted as contended.
- 5-second window is the rough hash-window assumption (a typical batch of UI clicks settling into one analyze cycle).
- 2% threshold is calibrated for "users actually contend, not just edit the same exposure on different days."

Run this query manually when the sustained-use minimum (≥4 weeks AND ≥500 delta events across ≥10 exposures) is met. The output of this query is what gates R5b — paste the three numbers into the PR description before any R5.2 code lands.

If R5b ships, the work is described in §R5b below.

**Files:**
- Modify: `packages/HimalayaUI/src/routes_peaks.jl`
- Modify: `packages/HimalayaUI/src/routes_analysis.jl`
- Modify: `packages/HimalayaUI/src/routes_exposures.jl`
- Modify: corresponding test files

For each delta-shaped route (see spec table):

- [ ] **Step (per route): extract and compare `If-Match`**

  ```julia
  function require_inputs_hash_match(db, req, exposure_id)::Union{HTTP.Response, Nothing}
      if_match = HTTP.header(req, "If-Match", "")
      isempty(if_match) && return nothing  # backwards-compat window: accept missing
      stored = first(Tables.rowtable(DBInterface.execute(db,
          "SELECT analysis_inputs_hash FROM exposures WHERE id = ?", [exposure_id]))).analysis_inputs_hash
      stored = ismissing(stored) ? "" : String(stored)
      if String(if_match) != stored
          current_state = build_current_state(db, exposure_id)
          return HTTP.Response(409,
              ["Content-Type" => "application/json"],
              JSON3.write(Dict(:current_inputs_hash => stored,
                                :current_state => current_state)))
      end
      nothing
  end
  ```

  Use it as a guard at the start of each delta-shaped route handler.

- [ ] **Step: deprecation warning when If-Match missing**

  Log `@warn "delta-shaped route called without If-Match header" route=...`. After one release window, escalate to 428 (Precondition Required).

### Task R5.3: Frontend SSE subscriber

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`
- Add: `packages/HimalayaUI/frontend/test/sse.test.ts`
- Add: `packages/HimalayaUI/frontend/e2e/multiplayer.spec.ts`

- [ ] **Step 1: SSE subscriber in App.tsx**

  ```typescript
  useEffect(() => {
    const es = new EventSource('/api/events');
    es.addEventListener('curation', (e) => {
      const event = JSON.parse(e.data);
      // Self-echo filter: skip events authored by this client. NOTE: if two
      // tabs share the same X-Username (same lab user with two browsers open),
      // each tab sees its own edits via the optimistic update path AND filters
      // the server echo. Edits from the *other* tab also get filtered, so the
      // second tab won't update until manual refetch. Acceptable in lab use;
      // if it becomes a problem, switch to a per-tab client id sent on requests.
      if (event.actor === username) return;
      const id = event.entity_id as number;
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
    });
    return () => es.close();
  }, [username, qc]);
  ```

- [ ] **Step 2: Two-context Playwright test**

  Headless: open two browser contexts pointing at the same dev server, log in as Alice and Bob (via different `X-Username`), edit the same exposure's peaks, assert both succeed and both see each other's edits within ~1s. Under R5a's LWW model: if Alice and Bob edit the same field simultaneously, the later write wins (no 409, no retry — that's R5b's job).

- [ ] **Step 3: Remove `autoReanalyze` chain**

  The hooks at [queries.ts:103-152](../../packages/HimalayaUI/frontend/src/queries.ts) currently make every peak edit a two-round-trip (peak op + reanalyze). With server-side memoization (R3), a single mutation is enough — the server reanalyzes if needed (skipping findpeaks/indexpeaks via hash check when nothing changed) and broadcasts. Delete `autoReanalyze` and inline the relevant logic.

- [ ] **Step 4: Run all tests, commit**

  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
  cd packages/HimalayaUI/frontend && npm test && npm run e2e && npm run build && cd -
  git add -p
  git commit -m "feat: SSE-driven multiplayer with LWW conflict resolution

  Other users' edits propagate within ~1s via SSE invalidation.
  Conflicts are resolved last-writer-wins. If R4 instrumentation later shows
  contention is frequent, R5b adds optimistic-concurrency If-Match retry."
  ```

### Task R5b (deferred): If-Match + 409 retry on delta-shaped routes

**Gating:** runs only when R5.2's pre-committed threshold trips — ≥2% of delta-shaped curation events have another delta-shaped event from a different actor within the 5-second hash window preceding them. See R5.2 above for the full query and rationale. The whole point of R4's instrumentation is to make this decision data-driven; the threshold is fixed here so the decision can't be reverse-engineered after the work has started.

If R5b ships, the work is:
- Backend: `If-Match` check on delta-shaped routes (per spec table); 409 + current-state body on mismatch.
- Frontend: `useMutationWithRetry` wrapper for one-retry-on-409 with refetch between attempts.
- Tests: two-context Playwright test where Alice and Bob both add peaks to the same exposure; both succeed via retry path.
- Backwards compat: missing `If-Match` accepted with deprecation warning for one release window, then escalated to 428.

---

## Final verification

- [ ] **Cross-DB-state validation**

  Test the full `open_db` → `serve` → `analyze` → `serve` flow on each of:
  - Fresh DB (created by current schema)
  - Pre-R0 DB snapshot (production DB before this plan)
  - Partially-migrated DB (R0+R1 only; R0+R1+R2 only; etc.)

- [ ] **Pipeline test suite**

  ```bash
  julia --project=. -e 'using Pkg; Pkg.test()'  # core
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'  # backend
  cd packages/HimalayaUI/frontend && npm test && npm run e2e && npm run build  # frontend
  ```

  All test counts in CLAUDE.md ("Current state" section) should match or exceed: ≥452 Julia (HimalayaUI), ≥90 Julia (core), ≥174 Vitest, ≥16 Playwright E2E. R5 adds ≥1 multiplayer Playwright test.

- [ ] **Smoke test the actual app**

  ```bash
  bin/himalaya serve /path/to/test/experiment --port 8080
  ```

  Open two browser tabs, edit peaks in one, observe the other update via SSE within ~1s. Toggle exclude on a peak, observe the StaleIndicesBanner condition resolve correctly via hash comparison.

- [ ] **Update CLAUDE.md "Current state"**

  Reflect new test counts and add a brief description of multiplayer + instrumentation:
  ```markdown
  - HimalayaUI Plan 7 — Multiplayer + Instrumentation Foundation:
    - Diff-based reanalysis preserves auto peak IDs.
    - Curation lives in dedicated `peak_curations` table, separate from machine output (`auto_peaks`).
    - Content-hash memoization skips findpeaks/indexpeaks when inputs unchanged.
    - Curation events recorded with structured payloads in `user_actions` (promoted to source-of-truth log).
    - SSE-driven multiplayer with If-Match conflict resolution on delta-shaped routes.
  ```

- [ ] **Commit final docs update**

  ```bash
  git add CLAUDE.md
  git commit -m "docs: bump current-state for multiplayer + instrumentation"
  ```

---

## Open items to revisit during implementation

- **Hash function: SHA-256 vs. cheap 64-bit.** This plan uses SHA-256 for stability. If profiling reveals trace-hash compute is hot (multi-megabyte files × thousands of exposures during reingest), reconsider for trace specifically.
- **Whether to ship R5b at all.** Gated on R4 instrumentation. If `analyze_run` events show concurrent edits to the same exposure are rare in practice, LWW under R5a is sufficient and R5b is wasted effort. Decision is downstream of data the system doesn't yet collect.
- **Oxygen.jl streaming API.** Task R5.0 spike resolves the concrete API choice before R5.1 commits. If Oxygen can't compose with SSE cleanly, fallback is a separate `HTTP.serve` for the event channel — the spec/plan stay valid; only the wiring code changes.
- **Speculative index UX under multiplayer.** When Alice's peak edit causes Bob's speculative index to rebind to a different peak, Bob's UI silently reflects the new state. A "your speculative was modified by Alice's edit" toast would help — out of scope here; track as a follow-up issue.
- **Reconnect backfill via `Last-Event-ID`.** Defer to R5a v2; clients re-sync via TanStack Query refetch on reconnect for now. The same mechanism would also cover the crash-during-broadcast window.
- **Authentication.** `X-Username` remains trusted. Multiplayer raises the impact of spoofing; flagged in spec.
- **R2 vs. `created_by`-only fallback.** Spec Open Question 7 covers this. If R2 implementation reveals unexpected complexity, the fallback (keep `peaks` table, add `created_by`) is a real escape hatch. R1+SSE+`created_by` ships ~80% of the multiplayer story without R2's table split, at the cost of leaving the snapshot/restore dance in `_persist_analysis_inner!` and complicating R4 payloads.

## What this plan does NOT do

- Does not migrate `indices`, `index_peaks`, `index_groups`, `index_group_members` to a derived/computed model.
- Does not implement any analysis-feedback loops (Bayesian phase priors, low-confidence flags, etc.). The instrumentation infrastructure is built; the loops it enables are deferred.
- Does not add offline editing, time-travel UI, or replay-based regression testing UI. The capabilities are enabled by R4; the surfaces are deferred.
