# Ingestion Redesign — Production Migration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Upgrade a **pre-ingestion-rework production database** (old schema, manifest-era grouping, zero loads) in place so it opens and renders fully on the new app — experiments show loads, grouping, and geometry — **with every existing analysis and curation preserved**. Concrete success: the dev-db copy (3 experiments / 139 samples / 678 exposures) opens, every experiment surface renders real data, and every exposure-keyed curation row count is identical before and after.

**Architecture:** Two halves, run by one operator-invoked CLI command.
1. **Schema half** (automatic, in `open_db` → `migrate_schema!`): the existing sentinel-gated migrations already do the additive schema upgrade. This plan changes exactly one — `migrate_exposures_experiment_id!`, whose `(experiment_id, filename)` collision `error()` hard-fails on real data — replacing the `error()` with a safe delete-only dedup.
2. **Data half** (explicit, one-shot CLI): a new `regroup_experiment!` reuses the Phase B/C **pure derivation** (`scan_directory` → `group_into_samples` → `derive_geometry`) and makes the DB match the **auto-grouping partition exactly as a fresh scan would** — one sample per `(load, slot)` cell — while **reusing existing rows** so human names + notes + curation survive. It never deletes/re-inserts an exposure (curation rides the stable exposure id).

**Core principle — the migrated DB is structurally identical to a freshly-scanned one.** Same grouping rules, same partition. The migration's only difference from a fresh scan is that it *relinks existing exposures* (preserving their ids/curation) instead of inserting new ones, and *reuses existing sample rows* (preserving human names/notes) instead of minting all-new ones. A reshot specimen splits into per-load samples just as it would on a fresh scan; re-uniting reshoots is the deferred reshoot-grouping capability (`docs/future-feature-ideas.md`), applied uniformly to migrated and new data.

**Why not just rescan:** `scan_and_group!` is insert-only — it *skips* existing exposures (filename already present, `ingest.jl:158`) and mints a parallel set of empty loads + new samples, stranding the old analyzed rows at `load_id=NULL`. So the rescan *derivation* is reusable but its *persist* path is not; the migration needs a retrofit persist that adopts existing rows.

**Tech Stack:** Julia, SQLite.jl / DBInterface.jl, stdlib `Test`. Backend at `packages/HimalayaUI/`. CLI in `packages/HimalayaUI/src/cli.jl`.

**Anchors verified 2026-06-20** against branch `ingestion-redesign` via two agent workflows (recon + 5-agent spec validation). Line numbers drift — confirm each with a quick `grep`/read before editing.

---

## Settled design decisions (do not re-litigate)

From the planning discussion (2026-06-20); fixed inputs.

- **B — Carry sample metadata.** Notes, `sample_tags`, `sample_messages`, `series_samples` are preserved: for a cell that reuses an old row they stay put; for a displaced old sample they are carried to the cell that absorbed its exposures.
- **C — Derive new geometry.** Legacy bare `energy_kev`/`flight_path_m` (source `'default'`) are overwritten by PRP/setup-derived values (`_update_geometry_if_not_user!` skips only `*_source='user'`). **Flag (don't block)** large discrepancies — surface `derive_geometry`'s discrepancy vector + old-vs-new deltas in the summary.
- **D — CLI, one-time.** The data half is an explicit operator command, never auto-run in `open_db`.
- **E — Files present.** Assume `data_dir`/`analysis_dir` resolve to the original files in production (verified true for the dev-db on this host).
- **Reshoot representation — match the auto-grouping (split).** A reshot specimen's exposures fall in multiple `(load, slot)` cells, so it becomes **multiple samples**, one per cell — exactly what a fresh scan produces (`samples.load_id` is single-valued; the auto-grouping partition is authoritative). The cell holding the specimen's earliest exposures **reuses the old sample row** (keeps its human name + notes); the later cell(s) become **new auto-named samples**. We do **not** re-unite them — that is the deferred reshoot-grouping capability (`docs/future-feature-ideas.md` → "Reshoot grouping for new scans"), and it is the *same* `merged_into_id` action for migrated and new data.
- **Dedup is delete-only.** The 30 `(experiment_id, filename)` collisions all share one `image_path`; non-survivor FK children are proven-redundant duplicates of the survivor's. No FK repoint — delete children (correct order) then the row.

---

## Ground truth (verified against the real dev-db, `regress.db`)

| Fact | Value | Consequence |
|---|---|---|
| Collision groups | 30, all pairs, all share one `image_path` | dedup safe & delete-only |
| Survivor rule | analyzed-first (`trace_hash NOT NULL`), tiebreak lowest `id` | never drops the only analyzed copy; `ORDER BY (trace_hash IS NULL), id` |
| Non-survivor FK children | 24 rows: none · 4 rows: `auto` index_groups + auto peaks (dup of survivor) · 2 rows: `rejection_reason=Flare` tag (survivor has it too) | delete children, no repoint |
| NULL-`sample_id` exposures | **0** | orphan `error()` (db.jl:554) won't fire; keep it as the fail-safe |
| Duplicate custom index_groups | **0** | `preflight_index_groups_uniqueness!` (db.jl:283) won't abort open_db |
| Sample partition vs `(load,slot)` | 108/136 single-cell (1:1 retrofit) · 28/136 multi-load (reshoots → split into per-load samples) | reconcile old partition → cell partition; reuse rows for names |
| Exposure nesting | `get_loads_rollup` nests samples under loads by `load_id`, exposures under samples by `sample_id` (`db.jl:2218-2231`) | each cell-sample renders under its own load; reshoot cells render faithfully per load |
| Source volumes | mounted, all 3 dirs reachable | derivation can run |

**Because the migrated partition matches a fresh scan,** the downstream gotchas of an anchor model do **not** arise: no sparse later-load (every cell is a real sample under its own load — verified: `_cluster_slots` returns ≥1 slot per non-empty load), no phantom cell on a later rescan (every cell already exists, so `scan_and_group!`'s `(experiment_id, load_index)` + `(experiment_id, load_id, slot_index)` dedup reuses it — verified against `ingest.jl:123,142`), and **the same `derive_sample_flags` output a fresh scan would produce** (the rollup→flag path is shared, so the migration introduces no flags beyond fresh-scan baseline). Note this last point is a *parity* guarantee, not a structural "zero flags" one — a cell with internal gaps in the ~0.5–1.5 mm window can still flag, exactly as on a fresh scan (real-data jitter ~0.3 mm makes it rare). The cost is that a reshot specimen appears as **two samples** until someone merges them — the deferred reshoot-grouping decision, applied consistently.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/db.jl` | `migrate_exposures_experiment_id!` dedup fix | MODIFY (Task 1) |
| `packages/HimalayaUI/src/ingest.jl` | new `regroup_experiment!` retrofit persist | MODIFY (Task 2) |
| `packages/HimalayaUI/src/cli.jl` | `upgrade-grouping` command (dry-run + apply) | MODIFY (Task 3) |
| `packages/HimalayaUI/test/test_migration_upgrade.jl` | new standalone test file | CREATE (all tasks) |
| `packages/HimalayaUI/test/runtests.jl` | test registry (GROUPS + ALL_ORDER) | MODIFY (Task 0) |
| `Makefile` (repo root) | `test-parallel` bucket list (line 38) | MODIFY (Task 0) |

**Out of scope:** the rest of the schema migrations (already correct), `scan_and_group!`'s insert path (untouched), any frontend change (the honest empty-states from `019602fc` cover the un-upgraded state), reshoot re-uniting (deferred).

---

## Conventions for every task

- **TDD.** Failing test → minimal implementation → verify → commit. One commit per task.
- **Run the single test file** during TDD: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_migration_upgrade.jl`
- **Legacy fixture, not `regress.db`.** Tests build a synthetic pre-rework DB in a tempdir (recipe in Task 0). Tests must not touch `/Volumes`.
- **Idempotent + transactional.** Every write is in `SQLite.transaction`; re-running is a no-op.
- **Commit trailers** (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01Uc6yH4xap9w3DcqRPF4v3V`.
- **Never `git add -A`** — stage explicit files only. Do not re-resolve `Manifest.toml` (SQLite pins deliberate).

---

## Task 0: Test harness + legacy-fixture builder

The hard part: `create_schema!` on this branch still emits the **legacy** `samples` shape (`name`+`display_name`+`samples_unique_name`) and omits the Phase-A additions; `migrate_schema!` is what upgrades. So a faithful pre-rework fixture is `SQLite.DB(path)` → `create_schema!` (**not** `open_db`, which also migrates) → hand-INSERT old-shape rows → close → reopen → `migrate_schema!`.

- [ ] Create `test/test_migration_upgrade.jl`. Register in **both** lists in `runtests.jl` (a load-time guard `error()`s the suite otherwise) and in the `Makefile`:

```julia
# runtests.jl — add a 'migration' bucket to GROUPS (every bucket carries test_http.jl for helper includes):
("migration", ["test_http.jl", "test_migration_upgrade.jl"]),
# AND add to ALL_ORDER:
"test_migration_upgrade.jl",
# Makefile:38 — append the bucket:
GROUPS := db pipeline routes events wire migration
```
(Guards: `runtests.jl:76-79` symdiff(GROUPS, ALL_ORDER) and `:85-90` on-disk-orphan both `error()` at load if you add to only one.)

- [ ] `build_legacy_db(dir)` — follow the Phase A `test_ingestion_schema.jl:274-313` recipe. Hand-write the old shape (`data_dir`/`analysis_dir` are NOT NULL — point them at fixture triplet dirs):

```julia
db = SQLite.DB(path); HimalayaUI.create_schema!(db)
# experiments(data_dir/analysis_dir → fixture dirs), samples WITH display_name,
# exposures with sample_id + bare-stem filename + NO experiment_id (old shape).
# Collision fixture: two exposures sharing (filename, image_path), one with trace_hash + an
#   auto index_groups/auto_peaks row, one empty. Reshoot fixture: one old sample whose stems
#   fall in two derived loads. Many-old→one-cell fixture: two old samples whose stems land in
#   one derived (load,slot) cell.
SQLite.close(db)   # then SQLite.DB(path) → migrate_schema! to exercise the migration
```

- [ ] For Task 1 dedup tests, reuse `_revert_exposures_denorm!` (`test_ingestion_schema.jl:103-112`): take a fully-migrated DB, un-apply the exposures migration (deletes the sentinel, restores the old `(sample_id,filename)` index, NULLs `experiment_id`), inject collisions, re-run `migrate_schema!`.
- [ ] On-disk fixture files (for Task 2 derivation) — reuse `test_ingestion_core.jl`'s `write_prp`/`write_setup_info`; **each stem needs a `.tif` (scan_directory's enumeration key), `.prp`, `.dat`, plus one `setup_info_*.txt` in analysis_dir.** The exposure's `filename` MUST equal the bare stem (no extension) or the retrofit's `WHERE filename = stem` finds nothing. Control `horizontal_position`/timestamps in the `.prp` to force the desired load/slot partition (so the test deterministically produces 1:1, reshoot, and many-old→one-cell shapes).
- [ ] `count_curation(db)` → NamedTuple of counts for the **exact** tables (note plurals): `auto_peaks, peak_curations, indices, index_peaks, index_group_members, index_groups, assignments, assignment_members, exposure_tags, exposure_sources, series_members, comparison_members`. Cross-task invariant: identical pre/post (minus the dedup'd non-survivors' redundant children).

## Task 1: Dedup fix in `migrate_exposures_experiment_id!`

**Anchor:** replace the **dupes-`error()` block at `db.jl:567-577`** with the dedup. It MUST run *before* the `DROP INDEX … / CREATE UNIQUE INDEX exposures_unique_filename ON exposures(experiment_id, filename)` at `db.jl:580-582` — else the index-create throws on the still-present collisions. Keep everything else (transaction `:532`, backfill `:547-551`, orphan `error()` `:554-559`, index swap `:580-586`, `_record_migration!` `:587`).

Current block being replaced:
```julia
dupes = Tables.rowtable(DBInterface.execute(db, """
    SELECT experiment_id, filename, COUNT(*) AS n
      FROM exposures WHERE filename IS NOT NULL
     GROUP BY experiment_id, filename HAVING n > 1"""))
isempty(dupes) || error("… manual merge required …")    # ← REPLACE THIS
```

- [ ] New logic: find collision groups; for each whose rows **all share one `image_path`** (`COUNT(DISTINCT image_path) = 1`), pick the survivor (`ORDER BY (trace_hash IS NULL), id LIMIT 1`) and delete the non-survivors.
- [ ] **Keep a hard `error()`** for any group with `COUNT(DISTINCT image_path) > 1` (genuinely distinct files; absent in real data).
- [ ] **FK-safe child deletes (load-bearing).** FK enforcement is most likely **ON** during this migration (earlier migrations leave `PRAGMA foreign_keys=ON` in their `finally` blocks, and the pragma can't change inside the open `SQLite.transaction`). Delete the non-survivors' children in strict child→parent order before the exposure rows, scoped to the non-survivor id set — mirror `routes_experiments.jl:464-496` (none of these FKs are `ON DELETE CASCADE` except `assignment_members.index_id`):
  1. `index_peaks`, `index_group_members` (via `indices JOIN exposures`) — grandchildren, FIRST
  2. `assignment_members`, `assignments` (PK is `exposure_id`), `index_groups`, `indices`, `auto_peaks`, `peak_curations`
  3. `exposure_sources` — **two columns**: `averaged_exposure_id` OR `source_exposure_id`
  4. `exposure_tags`
  5. then `DELETE FROM exposures WHERE id IN (<non-survivors>)`
- [ ] Tests: (a) same-path collision migrates clean — survivor + its curation intact, `count_curation` drops by exactly the non-survivors' redundant children, no orphaned FK rows, open succeeds; (b) re-`migrate_schema!` is a no-op (sentinel recorded); (c) a **differing**-`image_path` collision still `error()`s; (d) an **orphan** (NULL-`sample_id`) exposure still trips the orphan `error()` (real data has 0 — don't silently swallow).

## Task 2: `regroup_experiment!` — apply the auto-grouping partition, retrofit rows

New function in `ingest.jl`: `regroup_experiment!(db, experiment_id) -> NamedTuple` (summary: loads created, cells, samples retrofitted-in-place, samples newly created, samples displaced+deleted, exposures relinked, exposures-with-no-on-disk-file, geometry fields + **discrepancies**, reshoot/many-old-one-cell counts).

- [ ] **Derive (read-only, mirror `ingest.jl:68-93`):** resolve `data_dir`/`analysis_dir`/patterns from the experiments row; `metas = scan_directory(...)`; **`isempty(metas) && return (status=:empty, …)`** (mirror `ingest.jl:73`); build `prp_paths` and `setup_files` (`readdir(analysis_dir)` for `setup_info_*.txt`); `geo, disc = derive_geometry(prp_paths, setup_files)` — **capture `disc`** (decision C); `result = group_into_samples(metas)`.
- [ ] **Build the cell list + `byfile`.** Each cell = a `(load_index, slot_index)` from `GroupedLoad`/`GroupedSample` with its auto-name + stems. `byfile :: Dict{String,…}` keyed on **`ge.stem`** (NOT filename; `GroupedExposure` has `.stem`, `exposures.filename` stores the bare stem):
```julia
cells = NamedTuple[]   # (load_index, slot_index, name, name_source, grouping_source, stems)
byfile = Dict{String,@NamedTuple{...}}()
for gl in result.loads, gs in gl.samples
    push!(cells, (load_index=gl.load_index, slot_index=gs.slot_index,
                  name=gs.name, name_source=gs.name_source,
                  grouping_source=gs.grouping_source, stems=[ge.stem for ge in gs.exposures]))
    for ge in gs.exposures
        byfile[ge.stem] = (load_index=gl.load_index, slot_index=gs.slot_index,
            prp_path=ge.prp_path, timestamp=ge.timestamp, exposure_time=ge.exposure_time_s,
            horizontal_position=ge.horizontal_position_mm, scan_id=ge.scan_id, frame_no=ge.frame_no)
    end
end
```
- [ ] **In one `lock(_DB_WRITE_LOCK)` + `SQLite.transaction`:**
  1. **Create loads** — one row per derived load (dedup on `(experiment_id, load_index)`; reuse on re-run). **Set `frame_count = gl.frame_count`** (`db.jl:2192` → `ConfigurationBody.tsx:151` reads it; default 0 = empty timeline). Map `load_index => load_id`.
  2. **Assign each cell a sample row — greedy reuse, deterministic order** (cells sorted by `(load_index, slot_index)`). For each cell:
     - `candidates` = distinct `sample_id`s of the cell's existing exposures (`SELECT DISTINCT sample_id FROM exposures WHERE experiment_id=? AND filename IN (<cell stems>)`), excluding already-**claimed** ids.
     - If `candidates` non-empty: `owner` = lowest-id candidate; mark claimed; `UPDATE samples SET load_id=?, slot_index=?, name_source='user' WHERE id=owner` — **keep `name`** (reusing the old human label). (This is the 1:1 retrofit *and* a reshoot's earliest cell.)
     - Else: `owner = create_sample!(db; experiment_id, name=cell.name, grouping_source=cell.grouping_source, load_id, slot_index)` (already defaults `name_source='auto'` — no extra UPDATE needed). (A reshoot's later cell, or a cell whose old samples were all already claimed.) Duplicate names across cells are **allowed**: `samples_unique_name` is dropped by `migrate_samples_name_collapse!` (`db.jl:635`, "NO new UNIQUE on the label") — identity is `(load_id, slot_index)`, not the label — so no disambiguation is required for correctness.
     - Record `cell → owner`.
  3. **Relink exposures** — for each existing exposure with a `byfile` entry: `UPDATE exposures SET sample_id=<cell owner>, load_id=<cell load>, prp_path=?, timestamp=?, exposure_time=?, horizontal_position=?, scan_id=?, frame_no=? WHERE experiment_id=? AND filename=?`. **Marshal `Missing`→`nothing` and format timestamps** exactly like `ingest.jl:167-172` (`ismissing(x) ? nothing : x`; `Dates.format(ts,"yyyy-mm-ddTHH:MM:SS")`). Exposures with no `byfile` entry (file gone): leave, count.
  4. **Displaced old samples** — an old sample whose every cell was claimed by a lower-id sibling now holds **no** exposures. **Capture the `old sample_id → absorbing owner` map during step 3** (while relinking, before `sample_id` is overwritten — it can't be reconstructed afterward; tally which owner received the most of each old sample's exposures). Then detect: `SELECT id FROM samples WHERE experiment_id=? AND merged_into_id IS NULL AND id NOT IN (SELECT DISTINCT sample_id FROM exposures WHERE sample_id IS NOT NULL)`. For each, carry metadata to its absorbing owner — `series_samples`/`sample_tags`/`sample_messages` via the merge re-point UPDATE (`routes_grouping.jl:191-203`), and **`notes` via a separate same-row UPDATE** (`UPDATE samples SET notes = COALESCE(notes, ?) WHERE id=<owner>`) because that block does NOT touch `notes` (a column on the samples row). Then `DELETE` the empty row. Count + report. (Reshoots are NOT displaced — their old row was reused for the earliest cell; only genuinely-superseded many-old→one-cell rows are. In the real dev-db this count is ≈0.)
  5. **Geometry** — `_update_geometry_if_not_user!(db, experiment_id, geo)` (skips `*_source='user'`, overwrites `'default'`; runs inside a transaction).
  6. Set `experiments.ingest_status='complete'`, `last_scanned_at` = now. (Drop `scan_signature` — write-only/inert.)
- [ ] **No exposure is ever deleted/re-inserted** here (only `UPDATE`/relink) → exposure-keyed curation preserved. The only sample deletes are displaced-empty rows, whose metadata is carried first.
- [ ] Tests (Task 0 fixtures): 1:1 sample → retrofitted in place, same id, name kept, `load_id`/`slot_index` set; **reshoot** → TWO samples (earliest cell keeps the old name + id; later cell is a new `'auto'`-named row), each under its own load, `get_loads_rollup` shows each load with its cell's exposures; **many-old→one-cell** → lowest-id old sample owns the cell, the other is displaced (its notes/tags carried) and deleted; `count_curation` identical (exposure-keyed) before/after; every file-present exposure has non-NULL `load_id`+correct `sample_id`; geometry populated + discrepancies surfaced; **re-running is a no-op** (loads dedup, owners already claimed by id, UPDATEs idempotent); empty-metas returns `:empty`, writes nothing.

## Task 3: `upgrade-grouping` CLI command

There is **no `reingest` command** (manifest ingestion was deleted). Mirror `cli_serve` (`cli.jl:249-266`) — central DB, no positional path. Reuse `_resolve_experiment` (`cli.jl:17-52`) for `--experiment` (accepts nothing=sole/numeric-id/path/name). Add the branch to `main()` (`cli.jl:321-343`).

- [ ] `cli_upgrade_grouping(args)`: flags `--apply` (default **dry-run**) and `--experiment <id|name>` (default: all). `db = open_db(default_db_path())` (runs the schema half incl. Task 1 dedup). For each target experiment: if `isdir(data_dir)`, run `regroup_experiment!`; print a per-experiment table (dedup'd, loads + frame counts, samples retrofitted/created/displaced, exposures relinked, geometry fields + **discrepancies/deltas**, `count_curation` before/after). Unreachable dirs: skip + report.
- [ ] **Dry-run = rollback-by-exception** (no rollback flag exists; `SQLite.transaction` commits on normal return, rolls back on throw). For dry-run, run the work inside the transaction and `throw` a private `_DryRunRollback` sentinel at the end of the block, catch at the call site, print the computed summary. `--apply` lets it commit (one transaction per experiment → a failure isolates to one experiment).
- [ ] Print a leading "back up the DB file before `--apply`" reminder.
- [ ] Tests: dry-run writes nothing (DB byte-identical after); `--apply` produces the Task 2 end state; `--apply` twice is a no-op.

## Task 4: End-to-end verification

- [ ] **Round-trip test:** `build_legacy_db` (collisions + curation + reshoot + many-old→one-cell) → `migrate_schema!` (schema half + dedup) → `upgrade-grouping --apply` → assert: opens clean, `loads>0` with frame counts, partition matches the derivation (one sample per cell), all file-present exposures relinked, `count_curation` identical (exposure-keyed; minus the 30 dedup'd non-survivors' redundant children), geometry derived. **Do NOT assert** `q_units`/`q_units_source` populated (no writer — stays NULL/`'default'`), nor `*_source` non-`'default'` for fields the PRP/setup didn't supply.
- [ ] **Live render walk (operator runbook, in the test-file header):** copy the real dev-db, `upgrade-grouping --apply`, `serve`, walk every surface (experiment corpus/grouping/config + legacy samples/focus/series/loupe) confirming real data renders, 0 console errors. Verify-by-rendering gate; manual, not CI.
- [ ] **Idempotency check:** after `--apply`, run `cheap_change_check` (or hit the scan route) once and confirm `false`/no-op. Post-dedup `persisted` drops 30 but on-disk `.tif` count is unchanged (collisions share one file) → no-op, *unless* `data_dir` holds stray un-persisted `.tif` (calibration frames); if so investigate before declaring idempotent. (`serve()` arms no schedulers, so nothing auto-fires.) Confirm a manual rescan is a clean reuse (every cell already exists → no new/phantom samples).
- [ ] Record reshoot/displaced counts in the summary so we can eyeball them against the dev-db numbers (≈28 reshoots, ≈0 many-old→one-cell expected).

---

## Self-Review

- [ ] Dedup runs in place of the `:567-577` error block, before the index-create at `:580`; deletes only proven-redundant rows in child→parent order; survivor's curation untouched; differing-path + orphan still `error()`.
- [ ] No exposure deleted/re-inserted except the 30 non-survivors → curation invariant holds; only displaced-empty samples are deleted (metadata carried first).
- [ ] `regroup_experiment!` keys on `ge.stem`; sets `loads.frame_count`; reuses old rows greedily by lowest id (keeps human names); creates `'auto'` rows for extra/displaced cells; marshals `Missing`/timestamps; surfaces geometry discrepancies; empty-metas early-return; idempotent + transactional; dry-run writes nothing.
- [ ] Partition matches a fresh scan (one sample per cell) → reshoots render per-load, no sparse loads, no phantom on rescan, and the same `derive_sample_flags` output a fresh scan would give (no migration-introduced flags).
- [ ] Full backend suite green (`make test-parallel`) with the new `migration` bucket; verification asserts only writable columns.

---

## Open considerations (track, not blockers)

- **Reshoots appear as separate samples** until merged — the deferred reshoot-grouping capability, now uniform across migrated and new data. The migration records no link between a reshoot's cells; a future `merged_into_id` UX re-unites them on demand.
- **Repeated names are allowed.** `samples_unique_name` is dropped by `migrate_samples_name_collapse!` (`db.jl:635`, "NO new UNIQUE on the label"); identity is `(load_id, slot_index)`, so a new cell sample may share a name with a retained manifest sample without the INSERT throwing. An optional cosmetic coord-suffix could keep names visually distinct, but it is not required for correctness.
- **A reshoot's series recipe narrows.** A `series_samples` (recipe) entry pointing at a reshot specimen now resolves to the retained earliest-cell sample; its later-cell exposures (relinked to a new sample) drop out of that recipe entry. No FK breakage or data loss — the plate (`series_members`) is exposure-keyed and unaffected — and it is exactly what a fresh scan would have produced. The deferred reshoot-merge restores it.
- **`idempotent_responses`/`user_actions`:** the regroup writes no event/idempotent rows, so no purge is needed (unlike `migrate_samples_naming!`). The Task 1 dedup deletes exposures but `idempotent_responses` is not exposure-keyed → no orphaning.
