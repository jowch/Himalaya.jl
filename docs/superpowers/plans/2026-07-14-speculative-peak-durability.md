# Speculative-Index Peak Durability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make speculative ("custom") index peak assignments survive reanalysis, unify the stored basis convention, heal existing peak-less prod rows, and make the peak-less state legible in the UI.

**Architecture:** A new `speculative_peak_intents` table records the user's ratio→q assignments at creation and is never wiped; `_persist_analysis_inner!` re-resolves from intents instead of a join against live peak tables, and no longer commits destructive empty results. A one-shot sentinel migration rescales cubic speculative basis to the normalized convention, backfills intents, and NULLs `analysis_inputs_hash` on affected exposures so one reanalysis heals them. Frontend gets an "unresolved" chip and a `ci: 0` fix for the Miller regression band.

**Tech Stack:** Julia (SQLite.jl / DBInterface / Oxygen), React + TypeScript + Observable Plot, stdlib Test + Vitest.

**Spec:** `docs/superpowers/specs/2026-07-14-speculative-peak-durability-design.md` (read it first).

## Global Constraints

- The Julia backend suite is slow (5–10 min). Run **single test files** during TDD: from `packages/HimalayaUI/`, `julia --project=. -e 'using HimalayaUI; include("test/<file>.jl")'`. Full-suite runs: capture once to a per-branch log file and grep it (see `packages/HimalayaUI/test/AGENTS.md`).
- Vitest must run one-shot (`npx vitest run <file>` from `packages/HimalayaUI/frontend/`); a settings hook rejects watch mode.
- `npm run build` (from `packages/HimalayaUI/frontend/`) must pass before PR. Do **not** expect a green bare `npx tsc --noEmit` — `test/` has ~183 pre-existing errors; `npm run build` (src-only) is the gate.
- Stored basis convention (after this plan): `indices.basis` = q of the first ratio position for **every** kind (`predicted_q = basis × phaseratios(P; normalize=true)`).
- The `status = 'stale'` enum value is dead — `db.jl` normalizes it to `'candidate'` on every open. Never write it.
- New frontend colors come from existing `@theme` tokens (use `text-error`); stable selectors use `data-*` attributes, never Tailwind classes or copy text.
- All new SQL writes must run inside the existing transaction at each site (POST route's `with_idempotency` tx, `persist_analysis!`'s tx, the migration's own tx). Never open a nested transaction.

---

### Task 1: Normalize the basis convention at creation

**Files:**
- Modify: `packages/HimalayaUI/src/speculative.jl` (function `build_speculative_index`, ~lines 93-141)
- Test: `packages/HimalayaUI/test/test_speculative.jl` (append testset)

**Interfaces:**
- Consumes: `Himalaya.phaseratios(P; normalize=true)`, `Himalaya.Index{P}`, `Himalaya.fit`, `Himalaya.score` (all existing).
- Produces: `build_speculative_index(peak_rows, P, ratio_to_peak_id)` returns `(; basis, score, r_squared, lattice_d, residuals)` where `basis` is now the **normalized** basis. `insert_speculative_index!` (unchanged this task) stores it into `indices.basis`. Tasks 2, 3, 5 rely on this convention.

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_speculative.jl`:

```julia
@testset "speculative basis uses the normalized convention (cubic)" begin
    # Pn3m's first un-normalized ratio is √2. Before the fix the stored basis
    # was fit against un-normalized ratios, so predicted_q came back shrunk
    # by √2. The invariant: predicted_q[anchor_rpos] must equal the anchor q.
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")

    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.1)", [e_id])
    p1 = Int(DBInterface.lastrowid(res))

    # 1-peak (anchor-only) fixture: the LSQ fit through one point is exact,
    # so the assertion below is exact. A multi-peak fixture would deviate by
    # the fit residual — do not tighten this test with more peaks.
    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Pn3m,
        Dict{Int,Int}(1 => p1))

    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT basis, r_squared FROM indices WHERE id = ?", [new_id]))[1]
    predicted = HimalayaUI.predicted_q_for_phase("Pn3m", Float64(row.basis))
    @test predicted[1] ≈ 0.1 atol=1e-9   # pre-fix: 0.1/√2 ≈ 0.0707
    # 1-peak fit has zero residual DOF ⇒ R² is NaN, bound as NULL by SQLite.
    @test ismissing(row.r_squared)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_http.jl"); include("test/test_speculative.jl")' 2>&1 | tail -20
```
Expected: the new testset FAILS on `predicted[1] ≈ 0.1` (actual ≈ 0.0707). All pre-existing testsets in the file still pass (Lamellar's first ratio is 1, so its basis is identical under both conventions).

- [ ] **Step 3: Switch the fit to normalized ratios**

In `packages/HimalayaUI/src/speculative.jl`, `build_speculative_index` currently reads (after the sparse-vector construction):

```julia
    # Least-squares fit through assigned (ratio, q) pairs (intercept fixed at 0).
    # Mirrors `Himalaya.fit` exactly — extracts (idx, q) from the sparse vector.
    observed_ratios_full = Himalaya.phaseratios(P)  # un-normalized
    observed_ratios_used = observed_ratios_full[rpos_sorted]
    basis_unnorm = observed_ratios_used \ qvals  # 1/d in fit's terms

    idx = Himalaya.Index{P}(basis_unnorm, peaks_sv, sharpness_sv)
    fit_result = Himalaya.fit(idx)
    s          = Himalaya.score(idx)

    residuals = Dict{Int, Float64}()
    ratios_unnorm = observed_ratios_full
    for (rpos, qv) in zip(rpos_sorted, qvals)
        ideal = ratios_unnorm[rpos] * basis_unnorm
        residuals[rpos] = abs(qv - ideal)
    end

    (; basis = basis_unnorm,
       score = s,
       r_squared = fit_result.R²,
       lattice_d = fit_result.d,
       residuals = residuals)
```

Replace with (note: `ratios` — the normalized list — is already defined at the top of the function as `ratios = Himalaya.phaseratios(P; normalize = true)`):

```julia
    # Least-squares fit through assigned (ratio, q) pairs (intercept fixed at
    # 0), against NORMALIZED ratios so the stored basis means "q of the first
    # ratio position" — the same convention auto indices use (core:
    # predictpeaks = basis × phaseratios(normalize=true)). Every consumer
    # (predicted_q_for_phase, MillerPlot) assumes this scale; fitting against
    # un-normalized ratios shrank cubic predictions by the first raw ratio.
    observed_ratios_used = ratios[rpos_sorted]
    basis = observed_ratios_used \ qvals

    # Index basis is only carried, never read, by fit/score (fit refits d
    # internally from peaks; score reads peaks + sharpness) — but pass the
    # normalized value so the constructed Index is convention-correct.
    idx = Himalaya.Index{P}(basis, peaks_sv, sharpness_sv)
    fit_result = Himalaya.fit(idx)
    s          = Himalaya.score(idx)

    # Residuals are numerically identical under either convention
    # (ratios_unnorm·basis_unnorm ≡ ratios_normed·basis_normed).
    residuals = Dict{Int, Float64}()
    for (rpos, qv) in zip(rpos_sorted, qvals)
        ideal = ratios[rpos] * basis
        residuals[rpos] = abs(qv - ideal)
    end

    (; basis = basis,
       score = s,
       r_squared = fit_result.R²,
       lattice_d = fit_result.d,
       residuals = residuals)
```

- [ ] **Step 4: Run test to verify it passes**

Run: same command as Step 2.
Expected: ALL testsets in `test_speculative.jl` PASS (the new one plus the pre-existing Lamellar ones — `basis ≈ 1.0` still holds because Lamellar's normalized and raw ratios coincide).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/speculative.jl packages/HimalayaUI/test/test_speculative.jl
git commit -m "fix: speculative basis stored in the normalized convention"
```

---

### Task 2: Normalize the re-attach fits in the pipeline

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (speculative re-attach block, ~lines 396-516)
- Test: `packages/HimalayaUI/test/test_speculative.jl` (append testset)

**Interfaces:**
- Consumes: `_spec_synthetic_exposure(qs)` and `_spec_run_reanalyze!(db, exposure_id)` — file-local test helpers already defined in `test_speculative.jl` (build a curation-peak-only exposure and drive `_persist_analysis_inner!` directly).
- Produces: the re-attach loop's `basis_for_snap` and `new_basis` are normalized-basis values; `UPDATE indices SET basis = ...` on the success path stores the normalized convention. Task 5 rewrites parts of this block and assumes these fits are already normalized.

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_speculative.jl`:

```julia
@testset "re-attach keeps the normalized basis convention (cubic)" begin
    # Pn3m normalized ratios: [1, √3/√2, √4/√2, ...] ≈ [1, 1.2247, 1.4142, ...].
    # Place peaks at basis 0.1 × those ratios; rp3 initially unassigned.
    r = Himalaya.phaseratios(Himalaya.Pn3m; normalize = true)
    fx = _spec_synthetic_exposure([0.1, 0.1 * r[2], 0.1 * r[3]])
    p1, p2 = fx.peak_ids[1], fx.peak_ids[2]

    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id, Himalaya.Pn3m,
        Dict{Int,Int}(1 => p1, 2 => p2))

    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    row = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis FROM indices WHERE id = ?", [spec_id]))[1]
    predicted = HimalayaUI.predicted_q_for_phase("Pn3m", Float64(row.basis))
    # Post-reanalyze the recomputed basis must still be normalized: the first
    # predicted position sits at the anchor q. Pre-fix, new_basis was fit
    # against un-normalized ratios ⇒ predicted[1] ≈ 0.1/√2.
    @test predicted[1] ≈ 0.1 rtol=1e-6
    # And the mixed-scale discovery bug is gone: with a correct basis, rp3's
    # predicted position (0.1·r[3]) finds the third peak.
    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test [Int(x.ratio_position) for x in ip] == [1, 2, 3]
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_http.jl"); include("test/test_speculative.jl")' 2>&1 | tail -20
```
Expected: new testset FAILS on `predicted[1] ≈ 0.1`. (Discovery may or may not find rp3 pre-fix; the basis assertion is the reliable failure.)

- [ ] **Step 3: Switch all four re-attach fit/residual sites to normalized ratios**

In `packages/HimalayaUI/src/pipeline.jl`, inside the speculative re-attach loop:

3a. The per-index ratio setup currently reads:
```julia
            ratios_unnorm = Himalaya.phaseratios(P)
            ratios_normed = Himalaya.phaseratios(P; normalize = true)
            n             = length(ratios_normed)
```
Replace with (the un-normalized list has no remaining consumer after 3b-3d):
```julia
            ratios_normed = Himalaya.phaseratios(P; normalize = true)
            n             = length(ratios_normed)
```

3b. In the `basis_for_snap` block, both fit expressions
```julia
                ratios_unnorm[rpos_seed] \ qvals_seed
```
(one in the `length(ratio_to_peak) >= 2` branch, one in the `valid_snaps` branch) become:
```julia
                ratios_normed[rpos_seed] \ qvals_seed
```

3c. The success-path recompute currently reads:
```julia
            # Recompute basis/r²/d/score using new peak values
            observed_ratios_used = ratios_unnorm[rpos_sorted]
            new_basis = observed_ratios_used \ qvals
```
Replace with:
```julia
            # Recompute basis/r²/d/score using new peak values. Normalized
            # ratios — indices.basis means "q of the first ratio position"
            # for every kind (same convention as insert + auto indices).
            observed_ratios_used = ratios_normed[rpos_sorted]
            new_basis = observed_ratios_used \ qvals
```

3d. The residual write inside the re-insert loop currently reads:
```julia
                ideal = ratios_unnorm[rpos] * new_basis
```
Replace with:
```julia
                ideal = ratios_normed[rpos] * new_basis
```

After 3a-3d, `grep -n "ratios_unnorm" packages/HimalayaUI/src/pipeline.jl` must return nothing.

- [ ] **Step 4: Run test to verify it passes**

Run: same command as Step 2.
Expected: ALL testsets PASS — including the pre-existing "auto-discovers new peaks" and "revives from stale" Lamellar testsets (unchanged behavior when the first ratio is 1).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_speculative.jl
git commit -m "fix: pipeline speculative re-attach fits against normalized ratios"
```

---

### Task 3: Durability migration (table + basis rescale + backfill + heal trigger)

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (new constant + new migration function; call it at the END of `migrate_schema!`, after `migrate_experiment_config_label_to_name!(db)`)
- Create: `packages/HimalayaUI/test/test_migrate_speculative_durability.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` (add the include, alphabetical with the other `test_migrate_*` entries)

**Interfaces:**
- Consumes: `schema_migrations` sentinel pattern (`migrate_comparisons_to_series!` is the template, `db.jl:956+`), `resolve_phase` (`speculative.jl` — same module, resolved at call time), `Himalaya.phaseratios`, `comparison_now_iso()` (`db.jl`).
- Produces: table `speculative_peak_intents(index_id, ratio_position, q)` (PK `(index_id, ratio_position)`, FK → `indices(id) ON DELETE CASCADE`); sentinel name constant `MIGRATION_SPECULATIVE_PEAK_DURABILITY = "speculative_peak_durability"`. Tasks 4 and 5 read/write this table.

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/test/test_migrate_speculative_durability.jl`:

```julia
using Test, SQLite, DBInterface, Tables
using Himalaya

# Seed helper: one experiment/sample/exposure and one speculative index row
# written RAW (bypassing insert_speculative_index!) so we can simulate
# legacy-convention rows.
function _mig_fixture()
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")
    (; db, experiment_id = exp_id, exposure_id = e_id)
end

# Re-arm the one-shot migration so it can be driven directly against
# hand-seeded legacy rows (open_db already ran it on the fresh DB).
function _rearm!(db)
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
        [HimalayaUI.MIGRATION_SPECULATIVE_PEAK_DURABILITY])
end

@testset "migration: cubic speculative basis rescaled once" begin
    fx = _mig_fixture()
    # Legacy Pn3m row: un-normalized basis 0.1/√2 for an anchor at q=0.1.
    old_basis = 0.1 / first(Himalaya.phaseratios(Himalaya.Pn3m))
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Pn3m', ?, 'candidate', 'speculative')""",
        [fx.exposure_id, old_basis])
    ix_id = Int(DBInterface.lastrowid(res))

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    b = Float64(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis FROM indices WHERE id = ?", [ix_id]))[1].basis)
    @test b ≈ 0.1 atol=1e-12

    # Sentinel: a second direct call must be a no-op (no double rescale).
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)
    b2 = Float64(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis FROM indices WHERE id = ?", [ix_id]))[1].basis)
    @test b2 ≈ 0.1 atol=1e-12
end

@testset "migration: intents backfilled from surviving index_peaks" begin
    fx = _mig_fixture()
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)",
        [fx.exposure_id])
    p1 = Int(DBInterface.lastrowid(res))
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id,
        Himalaya.Lamellar, Dict{Int,Int}(1 => p1))
    # Simulate a pre-fix DB: no intents rows yet (Task 4 writes them at
    # creation; delete to model the legacy state).
    DBInterface.execute(fx.db,
        "DELETE FROM speculative_peak_intents WHERE index_id = ?", [spec_id])

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    intents = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position, q FROM speculative_peak_intents WHERE index_id = ?",
        [spec_id]))
    @test length(intents) == 1
    @test Int(intents[1].ratio_position) == 1
    @test Float64(intents[1].q) ≈ 0.05 atol=1e-12
end

@testset "migration: analysis_inputs_hash NULLed only for peak-less-speculative exposures" begin
    fx = _mig_fixture()
    # Exposure A: peak-less speculative (prod shape) → hash must be NULLed.
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Square', 0.17, 'candidate', 'speculative')""", [fx.exposure_id])
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'aaaa' WHERE id = ?",
        [fx.exposure_id])

    # Exposure B (same DB): healthy speculative → hash untouched.
    s2 = HimalayaUI.create_sample!(fx.db; experiment_id=fx.experiment_id, name="D2")
    e2 = HimalayaUI.create_exposure!(fx.db; sample_id=s2, filename="y")
    res = DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)", [e2])
    p1 = Int(DBInterface.lastrowid(res))
    HimalayaUI.insert_speculative_index!(fx.db, e2, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1))
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'bbbb' WHERE id = ?", [e2])

    _rearm!(fx.db)
    HimalayaUI.migrate_speculative_peak_durability!(fx.db)

    ha = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT analysis_inputs_hash AS h FROM exposures WHERE id = ?", [fx.exposure_id]))[1].h
    hb = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT analysis_inputs_hash AS h FROM exposures WHERE id = ?", [e2]))[1].h
    @test ismissing(ha)          # gate reopened
    @test String(hb) == "bbbb"   # healthy exposure untouched
end
```

Note: the backfill testset calls `insert_speculative_index!` then deletes intents — that DELETE is valid even before Task 4 exists (0 rows). The testsets are written to pass once Task 3 lands and stay green after Task 4.

- [ ] **Step 2: Register the test file and run to verify it fails**

Add to `packages/HimalayaUI/test/runtests.jl`, alphabetically among the existing includes:

```julia
include("test_migrate_speculative_durability.jl")
```

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_migrate_speculative_durability.jl")' 2>&1 | tail -20
```
Expected: FAIL — `UndefVarError: MIGRATION_SPECULATIVE_PEAK_DURABILITY not defined` (or `no such table: speculative_peak_intents`).

- [ ] **Step 3: Implement the migration**

In `packages/HimalayaUI/src/db.jl`, next to `MIGRATION_COMPARISONS_TO_SERIES` (top of file):

```julia
const MIGRATION_SPECULATIVE_PEAK_DURABILITY = "speculative_peak_durability"
```

Add the function (place it near `migrate_comparisons_to_series!`):

```julia
"""
    migrate_speculative_peak_durability!(db)

One-shot durability migration for speculative indices (spec:
2026-07-14-speculative-peak-durability-design.md):

1. Creates `speculative_peak_intents` — the durable record of the user's
   ratio→q assignments, written at creation and never wiped by analysis.
   Deliberately NOT in the `SCHEMA` const: `create_schema!` runs before
   `migrate_pk_to_autoincrement!` rename-rebuilds `indices`, and SQLite's
   RENAME tracking would rewrite this table's FK to `_migrate_old_indices`
   (same hazard `migrate_compare!` documents).
2. Rescales legacy speculative `basis` values from the un-normalized fit
   convention to the normalized one (factor = first raw phase ratio: √2
   Pn3m/Im3m, √3 Fm3m/Fd3m, √6 Ia3d; 1 elsewhere — no-op).
3. Backfills intents from surviving `index_peaks` rows (best available proxy
   for the creation-time assignment). Rows whose peak was deleted join to
   NULL q and are skipped — orphaned before this table existed, unrecoverable.
4. NULLs `exposures.analysis_inputs_hash` where a speculative index has zero
   `index_peaks` rows, reopening the `indexpeaks_skipped` memoization gate so
   the next analyze runs the re-attach loop (the heal pass). One-shot repair,
   not a standing gate bypass.

Steps 2-4 + the sentinel marker run in one transaction, sentinel-gated.
Must be called at the END of `migrate_schema!`: after
`migrate_pk_to_autoincrement!` + FK-heal (indices table settled), after
`migrate_r2_split_peaks!` (backfill join needs the post-R2 peak tables), and
after `migrate_series!` (creates `schema_migrations`).
"""
function migrate_speculative_peak_durability!(db::SQLite.DB)
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS speculative_peak_intents (
            index_id       INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
            ratio_position INTEGER NOT NULL,
            q              REAL    NOT NULL,
            PRIMARY KEY (index_id, ratio_position)
        )""")

    already = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM schema_migrations WHERE name = ?",
        [MIGRATION_SPECULATIVE_PEAK_DURABILITY]))
    isempty(already) || return nothing

    SQLite.transaction(db) do
        spec_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, phase, basis FROM indices WHERE kind = 'speculative'"))
        for r in spec_rows
            ismissing(r.basis) && continue
            P = resolve_phase(String(r.phase))
            P === nothing && continue
            factor = Float64(first(Himalaya.phaseratios(P)))
            factor == 1.0 && continue
            DBInterface.execute(db,
                "UPDATE indices SET basis = ? WHERE id = ?",
                [Float64(r.basis) * factor, Int(r.id)])
        end

        DBInterface.execute(db, """
            INSERT OR IGNORE INTO speculative_peak_intents (index_id, ratio_position, q)
            SELECT ip.index_id, ip.ratio_position, COALESCE(ap.q, pc.q)
            FROM index_peaks ip
            JOIN indices i ON i.id = ip.index_id AND i.kind = 'speculative'
            LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
            LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
            WHERE COALESCE(ap.q, pc.q) IS NOT NULL""")

        DBInterface.execute(db, """
            UPDATE exposures SET analysis_inputs_hash = NULL
             WHERE id IN (SELECT DISTINCT i.exposure_id FROM indices i
                           WHERE i.kind = 'speculative'
                             AND NOT EXISTS (SELECT 1 FROM index_peaks ip
                                              WHERE ip.index_id = i.id))""")

        DBInterface.execute(db,
            "INSERT INTO schema_migrations (name, applied_at) VALUES (?, ?)",
            [MIGRATION_SPECULATIVE_PEAK_DURABILITY, comparison_now_iso()])
    end
    nothing
end
```

Call it at the very end of `migrate_schema!`, immediately after `migrate_experiment_config_label_to_name!(db)`:

```julia
    # Speculative peak durability (2026-07-14 spec): intents table + legacy
    # basis rescale + heal trigger. Must run last — see the function docstring
    # for the three ordering constraints.
    migrate_speculative_peak_durability!(db)
```

- [ ] **Step 4: Run test to verify it passes**

Run: same command as Step 2. Expected: all three testsets PASS.
Also run `julia --project=. -e 'using HimalayaUI; include("test/test_db.jl")' 2>&1 | tail -5` — expected PASS (no interference with existing migrations).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/db.jl packages/HimalayaUI/test/test_migrate_speculative_durability.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat: speculative_peak_intents table + rescale/backfill/heal migration"
```

---

### Task 4: Write intents at creation

**Files:**
- Modify: `packages/HimalayaUI/src/speculative.jl` (function `insert_speculative_index!`, ~lines 173-218)
- Test: `packages/HimalayaUI/test/test_speculative.jl` (append testset)

**Interfaces:**
- Consumes: `speculative_peak_intents` table (Task 3).
- Produces: every speculative created through `insert_speculative_index!` has one intents row per assignment, `q` = the assigned peak's q at creation. Task 5's re-attach reads these.

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_speculative.jl`:

```julia
@testset "creation writes durable intents; index delete cascades them" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="x")
    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 1.0)", [e_id])
    p1 = Int(DBInterface.lastrowid(res))
    res = DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 2.0)", [e_id])
    p2 = Int(DBInterface.lastrowid(res))

    spec_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    intents = Tables.rowtable(DBInterface.execute(db,
        """SELECT ratio_position, q FROM speculative_peak_intents
           WHERE index_id = ? ORDER BY ratio_position""", [spec_id]))
    @test [(Int(r.ratio_position), Float64(r.q)) for r in intents] == [(1, 1.0), (2, 2.0)]

    # ON DELETE CASCADE (FK enforcement is ON via open_db).
    DBInterface.execute(db, "DELETE FROM indices WHERE id = ?", [spec_id])
    @test isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM speculative_peak_intents WHERE index_id = ?", [spec_id])))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_http.jl"); include("test/test_speculative.jl")' 2>&1 | tail -20
```
Expected: new testset FAILS (0 intents rows). Note: Task 3's backfill testset deletes creation-written intents before re-arming, so it stays green after this task.

- [ ] **Step 3: Write intents in `insert_speculative_index!`**

In `packages/HimalayaUI/src/speculative.jl`, the insert loop currently reads:

```julia
    for (rpos, peak_id) in ratio_to_peak_id
        pk_kind = _kind_for(db, exposure_id, peak_id)
        DBInterface.execute(db,
            """INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
               VALUES (?, ?, ?, ?, ?)""",
            [new_id, peak_id, pk_kind, rpos, built.residuals[rpos]])
    end
    new_id
```

Replace with (all writes stay after `build_speculative_index` — which validates every peak id and throws before any INSERT — inside the caller's `with_idempotency` transaction):

```julia
    q_by_id = Dict{Int, Float64}(Int(pr.id) => Float64(pr.q) for pr in peak_rows)
    for (rpos, peak_id) in ratio_to_peak_id
        pk_kind = _kind_for(db, exposure_id, peak_id)
        DBInterface.execute(db,
            """INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
               VALUES (?, ?, ?, ?, ?)""",
            [new_id, peak_id, pk_kind, rpos, built.residuals[rpos]])
        # Durable intent: index_peaks is the per-analysis resolved view (wiped
        # and rebuilt by _persist_analysis_inner!); this row is the user's
        # assignment and survives every wipe. Frozen at creation.
        DBInterface.execute(db,
            """INSERT INTO speculative_peak_intents (index_id, ratio_position, q)
               VALUES (?, ?, ?)""",
            [new_id, rpos, q_by_id[peak_id]])
    end
    new_id
```

- [ ] **Step 4: Run test to verify it passes**

Run: same command as Step 2. Expected: ALL testsets PASS.
Also run the migration file (its backfill fixture depends on this function):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_migrate_speculative_durability.jl")' 2>&1 | tail -5
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/speculative.jl packages/HimalayaUI/test/test_speculative.jl
git commit -m "feat: record speculative peak intents at creation"
```

---

### Task 5: Re-attach from intents, non-destructively

**Files:**
- Modify: `packages/HimalayaUI/src/pipeline.jl` (snapshot query ~line 273; re-attach loop ~lines 349-484)
- Modify: `packages/HimalayaUI/test/test_speculative.jl` (rewrite the "revives from stale" testset; append round-trip testset)
- Modify: `packages/HimalayaUI/src/AGENTS.md` (one bullet documenting the intents table)

**Interfaces:**
- Consumes: `speculative_peak_intents` rows (Tasks 3-4); `eff_lookup`/`eff_by_id`/`EXCLUDE_TOL` (existing in the function); normalized fits (Task 2).
- Produces: behavior contract for Task 6 and the frontend — a speculative with zero resolvable intents keeps `status='candidate'`, empty `index_peaks`, untouched score/basis, and intents intact; ≥1 resolved peak rebuilds rows and recomputes stats.

- [ ] **Step 1: Update the existing testset + write the new failing testset**

In `packages/HimalayaUI/test/test_speculative.jl`, replace the entire
`@testset "speculative revives from stale when matching peaks return"` block with:

```julia
@testset "speculative survives losing its peaks and re-attaches from intents" begin
    # Lamellar basis 0.05: peaks at 0.05 (rp1), 0.10 (rp2), 0.15 (rp3).
    fx = _spec_synthetic_exposure([0.05, 0.10, 0.15])
    p1, p2, p3 = fx.peak_ids[1], fx.peak_ids[2], fx.peak_ids[3]
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    # Delete both assigned peaks. Intents (rp1→0.05, rp2→0.10) survive; the
    # basis seed from intent q's lets discovery still snap p3 into rp3, so
    # the index stays live with 1 resolved peak (old behavior: stale + wiped).
    DBInterface.execute(fx.db,
        "DELETE FROM peak_curations WHERE id IN (?, ?)", [p1, p2])
    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    @test String(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status FROM indices WHERE id = ?", [spec_id]))[1].status) == "candidate"
    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position, peak_id FROM index_peaks WHERE index_id = ?", [spec_id]))
    @test length(ip) == 1
    @test Int(ip[1].ratio_position) == 3
    intents = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position FROM speculative_peak_intents WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test [Int(r.ratio_position) for r in intents] == [1, 2]  # intents untouched

    # Re-add peaks at the intent positions: next analysis re-resolves them.
    DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)",
        [fx.exposure_id])
    DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.10)",
        [fx.exposure_id])
    _spec_run_reanalyze!(fx.db, fx.exposure_id)
    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test [Int(r.ratio_position) for r in ip] == [1, 2, 3]
end
```

Append two new testsets:

```julia
@testset "zero-resolution analysis is non-destructive" begin
    # Only the two assigned peaks exist; delete both and give discovery
    # nothing to find (no other peaks) — 0 resolved.
    fx = _spec_synthetic_exposure([0.05, 0.10])
    p1, p2 = fx.peak_ids[1], fx.peak_ids[2]
    spec_id = HimalayaUI.insert_speculative_index!(fx.db, fx.exposure_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))
    pre = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT basis, score FROM indices WHERE id = ?", [spec_id]))[1]

    DBInterface.execute(fx.db,
        "DELETE FROM peak_curations WHERE id IN (?, ?)", [p1, p2])
    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    row = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status, basis, score FROM indices WHERE id = ?", [spec_id]))[1]
    @test String(row.status) == "candidate"       # never 'stale' (dead enum)
    @test Float64(row.basis) ≈ Float64(pre.basis) # last known stats kept
    @test Float64(row.score) ≈ Float64(pre.score)
    @test isempty(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT 1 FROM index_peaks WHERE index_id = ?", [spec_id])))
    @test length(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT 1 FROM speculative_peak_intents WHERE index_id = ?", [spec_id]))) == 2

    # Peaks return → heals on the next run, from intents alone.
    DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.05)",
        [fx.exposure_id])
    DBInterface.execute(fx.db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.10)",
        [fx.exposure_id])
    _spec_run_reanalyze!(fx.db, fx.exposure_id)
    @test length(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT 1 FROM index_peaks WHERE index_id = ?", [spec_id]))) == 2
end

@testset "prod-shape heal: peak-less speculative re-attaches from stored basis" begin
    # The 11 prod rows: basis set, no intents, no index_peaks. Discovery from
    # the stored (normalized) basis re-attaches everything it can.
    fx = _spec_synthetic_exposure([0.05, 0.10, 0.15])
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Lamellar', 0.05, 'candidate', 'speculative')""",
        [fx.exposure_id])
    spec_id = Int(DBInterface.lastrowid(res))

    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    row = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT status, score FROM indices WHERE id = ?", [spec_id]))[1]
    @test String(row.status) == "candidate"
    @test !ismissing(row.score)                   # stats recomputed
    ip = Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT ratio_position FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
        [spec_id]))
    @test [Int(r.ratio_position) for r in ip] == [1, 2, 3]
end
```

- [ ] **Step 2: Run to verify the new/changed testsets fail**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_http.jl"); include("test/test_speculative.jl")' 2>&1 | tail -30
```
Expected: the rewritten "survives losing its peaks" testset FAILS (current code marks stale + wipes; intents aren't read), "zero-resolution" FAILS (status flips are ephemeral but the wipe commits — `index_peaks` empty assertion passes, the intent-based heal assertion fails), "prod-shape heal" may already PASS (stored-basis fallback exists) — that one is a pin, not a change.

- [ ] **Step 3: Rewrite the snapshot source and the failure path**

In `packages/HimalayaUI/src/pipeline.jl`:

3a. Replace the snapshot query (currently `speculative_assignments = Tables.rowtable(DBInterface.execute(db, """SELECT i.id AS index_id, ip.ratio_position, ip.peak_kind, ip.peak_id AS old_peak_id, COALESCE(ap.q, pc.q) AS q_value FROM indices i JOIN index_peaks ip ... WHERE i.exposure_id = ? AND i.kind = 'speculative'"""))`) with:

```julia
    # The user's durable ratio→q assignments. index_peaks is the per-analysis
    # resolved view (wiped and rebuilt below); intents are never wiped, so a
    # failed re-resolution is retryable on every subsequent analysis.
    speculative_intents = Tables.rowtable(DBInterface.execute(db, """
        SELECT spi.index_id, spi.ratio_position, spi.q AS q_value
        FROM speculative_peak_intents spi
        JOIN indices i ON i.id = spi.index_id
        WHERE i.exposure_id = ? AND i.kind = 'speculative'
        """, [exposure_id]))
```

Also update the comment block above it (it references "Snapshot speculative indices' index_peaks rows by q-value") to say intents are the source.

3b. Inside `if !isempty(speculative_index_ids)`, replace the whole `_resolve_peak_id` closure with:

```julia
        # Match an intent q against the current effective peak set. Curation
        # peaks keep their exact q across analyses (exact-lookup hit); auto
        # peaks may drift slightly under re-detection (fuzzy fallback, same
        # tolerance the old snapshot resolution used).
        function _resolve_intent_peak(qv::Float64)
            info = get(eff_lookup, qv, nothing)
            info !== nothing && return info[1]
            tol_val = max(EXCLUDE_TOL, abs(qv) * 0.001)
            best, best_delta = nothing, Inf
            for i in eachindex(eff.q)
                d = abs(eff.q[i] - qv)
                if d < best_delta && d <= tol_val
                    best_delta = d
                    best = eff.peak_id[i]
                end
            end
            best
        end
```

3c. Replace the `by_index` grouping source:

```julia
        # Group intent rows by index_id.
        by_index = Dict{Int, Vector{NamedTuple}}()
        for r in speculative_intents
            push!(get!(by_index, Int(r.index_id), NamedTuple[]), r)
        end
```

3d. Replace the per-index snap resolution loop (currently `pid = _resolve_peak_id(s)`):

```julia
            snaps = get(by_index, ix_id, NamedTuple[])
            ratio_to_peak = Dict{Int, Tuple{Int, Float64, Float64}}()  # rpos → (peak_id, q, sharpness)
            for s in snaps
                rpos = Int(s.ratio_position)
                pid  = _resolve_intent_peak(Float64(s.q_value))
                pid === nothing && continue
                pr = get(eff_by_id, pid, nothing)
                pr === nothing && continue
                ratio_to_peak[rpos] = (pid, Float64(pr.q), Float64(pr.sharpness))
            end
```

3e. Simplify `basis_for_snap` (intent q is `NOT NULL`, so the old `ismissing(s.q_value)` filter is dead):

```julia
            basis_for_snap = if length(ratio_to_peak) >= 2
                # Nominal path: use currently-resolved peaks to refit basis.
                rpos_seed  = sort(collect(keys(ratio_to_peak)))
                qvals_seed = [ratio_to_peak[r][2] for r in rpos_seed]
                ratios_normed[rpos_seed] \ qvals_seed
            elseif length(snaps) >= 2
                # Recovery path: fit from the intent q's — they encode the
                # user's hypothesis even if the underlying peaks are gone.
                snap_pairs = sort([(Int(s.ratio_position), Float64(s.q_value)) for s in snaps])
                ratios_normed[[first(p) for p in snap_pairs]] \ [last(p) for p in snap_pairs]
            else
                # Last resort: the persisted (normalized) basis on the row.
                if ismissing(ix_row.basis)
                    nothing
                else
                    stored = Float64(ix_row.basis)
                    stored > 0 ? stored : nothing
                end
            end
```

3f. Replace the failure guard (currently `if length(ratio_to_peak) < 2 ... SET status = 'stale' ... continue`):

```julia
            # Non-destructive failure: zero resolved peaks leaves the resolved
            # view empty for this run — intents survive and the next analysis
            # retries from them. No status write: 'stale' is a dead enum value
            # (db.jl R3.2 normalizes it to 'candidate' on every open); the UI
            # signals this state from peaks.length === 0. One resolved peak is
            # enough to keep the index live — anchor-only speculatives are
            # legal at creation and stay legal here (their R² is NaN → NULL,
            # which the serializers already guard).
            isempty(ratio_to_peak) && continue
```

- [ ] **Step 4: Run tests to verify they pass**

Run: same command as Step 2.
Expected: ALL testsets in `test_speculative.jl` PASS — including the pre-existing real-data "pipeline preserves speculative across re-analyze" (intents resolve by exact/fuzzy q the same way the old snapshot did) and Task 2's cubic testset.
Then confirm no dangling references: `grep -n "_resolve_peak_id\|speculative_assignments" packages/HimalayaUI/src/pipeline.jl` → no matches.

- [ ] **Step 5: Document the invariant and commit**

Add one bullet to `packages/HimalayaUI/src/AGENTS.md` in whatever section covers the pipeline/indices (match surrounding style):

```markdown
- **`speculative_peak_intents` is the durable record; `index_peaks` is a resolved view.** Analysis wipes and rebuilds `index_peaks` for speculatives every run, resolving from intents (creation-time ratio→q assignments, frozen, cascade-deleted with the index). Never write intents from discovery results, and never let a failed re-resolution delete them — zero resolved peaks must leave the row `candidate` with empty `index_peaks` (do not write `status='stale'`; it's a dead enum normalized away on open).
```

```bash
git add packages/HimalayaUI/src/pipeline.jl packages/HimalayaUI/test/test_speculative.jl packages/HimalayaUI/src/AGENTS.md
git commit -m "feat: re-attach speculative peaks from durable intents, non-destructively"
```

---

### Task 6: SSE post_state contract test

**Files:**
- Test: `packages/HimalayaUI/test/test_speculative.jl` (append testset — test-only task)

**Interfaces:**
- Consumes: `HimalayaUI._serialized_indices_for_broadcast(db, exposure_id)` (existing, `pipeline.jl` — builds the `analyze_run` `post_state` indices payload), Task 5 heal behavior.
- Produces: a pin on the contract that healed speculative state reaches SSE clients through the existing plumbing (spec §Rollout; `docs/contract-testing.md`).

- [ ] **Step 1: Write the test**

Append to `packages/HimalayaUI/test/test_speculative.jl`:

```julia
@testset "analyze_run post_state carries healed speculative peaks (SSE contract)" begin
    # A peak-less speculative heals during analysis; the post_state payload
    # every SSE subscriber receives must contain its re-attached peaks —
    # this is the only channel multiplayer clients converge through.
    fx = _spec_synthetic_exposure([0.05, 0.10])
    res = DBInterface.execute(fx.db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, 'Lamellar', 0.05, 'candidate', 'speculative')""",
        [fx.exposure_id])
    spec_id = Int(DBInterface.lastrowid(res))

    _spec_run_reanalyze!(fx.db, fx.exposure_id)

    payload = HimalayaUI._serialized_indices_for_broadcast(fx.db, fx.exposure_id)
    spec = only(filter(d -> Int(d[:id]) == spec_id, payload))
    @test String(spec[:kind]) == "speculative"
    @test length(spec[:peaks]) == 2
    @test String(spec[:status]) == "candidate"
end
```

Note for the implementer: the payload Dicts are Symbol-keyed (built via `row_to_json`, same as the GET route). If the key access errors, read the tail of `_serialized_indices_for_broadcast` in `pipeline.jl` and match its key type — the assertion targets (kind/peaks/status) are the contract, the key spelling is not.

- [ ] **Step 2: Run the test**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using HimalayaUI; include("test/test_http.jl"); include("test/test_speculative.jl")' 2>&1 | tail -10
```
Expected: PASS immediately (Tasks 2+5 already made heal work; this pins the broadcast layer). If it fails on key spelling, fix the test's key access per the note — not the source.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_speculative.jl
git commit -m "test: pin analyze_run post_state carries healed speculative peaks"
```

---

### Task 7: "peaks unresolved" chip on IndexCard

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx` (`IndexCard`, ~lines 45-110)
- Test: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx` (append describe block)

**Interfaces:**
- Consumes: `IndexEntry.kind` / `IndexEntry.peaks` (existing API shape — unchanged by this plan).
- Produces: `data-unresolved` attribute on the card `<li>` and `data-testid="index-unresolved-<id>"` on the chip (stable selectors for Vitest/E2E).

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx` (uses the file's existing `mockAll` helper and `renderWithProviders`):

```tsx
describe("<PhasePanel> — unresolved speculative chip", () => {
  const specBase = {
    id: 30, exposure_id: 42, phase: "Hexagonal", basis: 0.19, score: null,
    r_squared: null, lattice_d: 38.2, status: "candidate", kind: "speculative",
    inputs_hash: null, predicted_q: [0.19, 0.33],
  };

  it("peak-less speculative shows the chip and data-unresolved", async () => {
    mockAll(
      [{ ...specBase, peaks: [] }],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const chip = await screen.findByTestId("index-unresolved-30");
    expect(chip).toHaveTextContent(/peaks unresolved/i);
    const li = document.querySelector('li[data-index-id="30"]');
    expect(li).toHaveAttribute("data-unresolved");
  });

  it("speculative with peaks keeps the count pill, no chip", async () => {
    mockAll(
      [{ ...specBase, peaks: [
        { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.19 },
      ] }],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await screen.findByText("Hexagonal");
    expect(screen.queryByTestId("index-unresolved-30")).toBeNull();
    expect(screen.getByText(/1 peaks/)).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npx vitest run test/PhasePanel.test.tsx 2>&1 | tail -15
```
Expected: the two new tests FAIL (`findByTestId` times out); existing tests PASS.

- [ ] **Step 3: Implement the chip**

In `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx`, `IndexCard`:

Add after the existing `const isSpeculative = index.kind === "speculative";`:

```tsx
  // Peak-less speculative: the durable intents exist server-side but no
  // current peak resolves them (spec: 2026-07-14-speculative-peak-durability).
  const unresolved = isSpeculative && index.peaks.length === 0;
```

Add `data-unresolved={unresolved || undefined}` to the `<li>`'s attributes (next to `data-speculative`).

Replace the peak-count pill:

```tsx
          <span className="ml-auto px-1.5 py-0.5 border border-border-soft rounded-full text-xs text-fg-dim shrink-0">
            {index.peaks.length} peaks
          </span>
```

with:

```tsx
          {unresolved ? (
            <span
              data-testid={`index-unresolved-${index.id}`}
              className="ml-auto px-1.5 py-0.5 border border-error rounded-full text-xs text-error shrink-0"
            >
              peaks unresolved
            </span>
          ) : (
            <span className="ml-auto px-1.5 py-0.5 border border-border-soft rounded-full text-xs text-fg-dim shrink-0">
              {index.peaks.length} peaks
            </span>
          )}
```

(`text-error` is an existing `@theme` token already used in this file — no new colors.)

- [ ] **Step 4: Run tests to verify they pass**

Run: same command as Step 2. Expected: ALL tests in the file PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx packages/HimalayaUI/frontend/test/PhasePanel.test.tsx
git commit -m "feat: unresolved-peaks chip on peak-less speculative index cards"
```

---

### Task 8: Suppress the NaN confidence band in MillerPlot

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MillerPlot.tsx` (regression mark options, ~lines 87-96)
- Test: `packages/HimalayaUI/frontend/test/MillerPlot.test.tsx` (append tests)

**Interfaces:**
- Consumes: `Plot.linearRegressionY` `ci` option (`@observablehq/plot@0.6.17`: "The confidence interval in (0, 1), or 0 to hide bands"; the regression line renders unconditionally).
- Produces: no `<path d="M…,NaN…">` SVG console errors for 2-point regressions.

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/frontend/test/MillerPlot.test.tsx` (the file already mocks `@observablehq/plot`; mirror its existing `linearRegressionY` assertion test):

```tsx
describe("<MillerPlot> — confidence band suppression", () => {
  const ix = (peaks: IndexEntry["peaks"]): IndexEntry => ({
    id: 5, exposure_id: 1, phase: "Hexagonal", basis: 0.1, score: 0.5,
    r_squared: 1.0, lattice_d: 38, ngc: null, status: "candidate",
    kind: "speculative", inputs_hash: null,
    predicted_q: [0.1, 0.173, 0.2],
    peaks,
  });

  it("passes ci: 0 for a 2-point regression (NaN band guard)", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.linearRegressionY as unknown as { mockClear: () => void }).mockClear();
    render(<MillerPlot indices={[ix([
      { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.1 },
      { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.173 },
    ])]} />);
    const calls = (Plot.linearRegressionY as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    expect(calls.length).toBeGreaterThan(0);
    expect((calls[0]![1] as { ci?: number }).ci).toBe(0);
  });

  it("leaves ci at its default for a 3-point regression", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.linearRegressionY as unknown as { mockClear: () => void }).mockClear();
    render(<MillerPlot indices={[ix([
      { peak_id: 1, ratio_position: 1, residual: 0, q_observed: 0.1 },
      { peak_id: 2, ratio_position: 2, residual: 0, q_observed: 0.173 },
      { peak_id: 3, ratio_position: 3, residual: 0, q_observed: 0.2 },
    ])]} />);
    const calls = (Plot.linearRegressionY as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    expect(calls.length).toBeGreaterThan(0);
    expect((calls[0]![1] as { ci?: number }).ci).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npx vitest run test/MillerPlot.test.tsx 2>&1 | tail -15
```
Expected: first new test FAILS (`ci` is `undefined`); second passes; existing tests pass.

- [ ] **Step 3: Implement**

In `packages/HimalayaUI/frontend/src/components/MillerPlot.tsx`, the regression mark push currently reads:

```tsx
      regressionMarks.push(
        Plot.linearRegressionY(rows, {
          x: "ratio",
          y: "q",
          stroke: color,
          strokeOpacity: opacity,
          strokeWidth: isHovered ? 1.5 : 1,
          ...(dashed ? { strokeDasharray: isSpeculative ? "2,3" : "4,3" } : {}),
        }),
      );
```

Replace with:

```tsx
      regressionMarks.push(
        Plot.linearRegressionY(rows, {
          x: "ratio",
          y: "q",
          stroke: color,
          strokeOpacity: opacity,
          strokeWidth: isHovered ? 1.5 : 1,
          // n=2 has zero residual DOF: the confidence band divides by (n-2)
          // and renders an all-NaN <path> (SVG console error on every draw).
          // ci: 0 hides only the band; the regression line always renders.
          ...(rows.length < 3 ? { ci: 0 } : {}),
          ...(dashed ? { strokeDasharray: isSpeculative ? "2,3" : "4,3" } : {}),
        }),
      );
```

- [ ] **Step 4: Run tests, then the frontend build gate**

```bash
npx vitest run test/MillerPlot.test.tsx 2>&1 | tail -5
npm run build 2>&1 | tail -5
```
Expected: tests PASS; build PASSES (this is the whole-plan frontend gate — run it here since this is the last frontend task).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MillerPlot.tsx packages/HimalayaUI/frontend/test/MillerPlot.test.tsx
git commit -m "fix: suppress NaN confidence band on 2-point Miller regressions"
```

---

### Task 9: Full-suite verification

**Files:** none (verification only).

- [ ] **Step 1: Run the full Julia backend suite, captured**

Run (from repo root; the suite is slow — capture once, grep the log):
```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test-speculative-durability.out 2>&1; tail -20 /tmp/jl-test-speculative-durability.out
```
Expected: all pass. Known intermittent flakes that are NOT regressions (re-run once if hit): fast-skip P99 latency, `GET /api/health` port race.

- [ ] **Step 2: Run the full Vitest suite, captured**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm test > /tmp/vitest-speculative-durability.out 2>&1; tail -10 /tmp/vitest-speculative-durability.out
```
Expected: all pass.

- [ ] **Step 3: Commit anything outstanding; done**

Working tree should already be clean. The prod heal (deploy → backup DB → `bin/himalaya analyze <dir>` per experiment → verify with the spec §D SQL) is an operator action after merge — it is in the spec runbook, not this plan.
