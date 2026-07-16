# packages/HimalayaUI/test/test_migration_upgrade.jl
#
# Task 0 of the ingestion-migration plan: the test harness + legacy-fixture
# builder consumed by later migration tasks (Task 1 dedup, Task 2
# regroup_experiment!, Task 4 round-trip). This file owns NO production-code
# behavior — it only provides reusable helpers + a smoke @testset proving the
# legacy fixture builds and migrate_schema! upgrades it cleanly.
#
# The hard part (see test_ingestion_schema.jl "legacy DB upgrades cleanly"):
# create_schema! on this branch still emits the LEGACY samples/exposures shape
# (samples WITH name+display_name; exposures with sample_id + bare-stem
# filename + an un-backfilled NULL experiment_id). migrate_schema! is what
# upgrades: it backfills experiment_id, collapses samples.display_name→name, and
# swaps the dedup key to (experiment_id, filename).
# So a faithful pre-rework fixture is:
#     SQLite.DB(path) -> create_schema! (NOT open_db, which also migrates)
#       -> hand-INSERT old-shape rows -> close
#       -> reopen -> migrate_schema!
#
# When run standalone (`julia --project=packages/HimalayaUI
# packages/HimalayaUI/test/test_migration_upgrade.jl`) the includes below pull
# in the shared fixture writers from test_ingestion_core.jl. Under the suite
# they are already in scope; the isdefined guards make the re-include a no-op.
#
# ===========================================================================
# OPERATOR RUNBOOK — the manual verify-by-rendering gate (Task 4, NOT CI)
# ===========================================================================
# These automated tests prove the migration on SYNTHETIC fixtures only (they
# never touch /Volumes). The capstone's "real data renders" guarantee is a
# manual gate an operator runs ONCE against a copy of the production dev-db
# before merging. Steps (do NOT automate — the source volumes are host-mounted):
#
#   1. Copy the real dev-db to a scratch path so the original is never mutated:
#        cp /path/to/prod/himalaya.db /tmp/upgrade-smoke.db
#   2. Dry-run first and eyeball the summary (NOTHING is committed):
#        HIMALAYA_DB_PATH=/tmp/upgrade-smoke.db bin/himalaya upgrade-grouping
#      Expect (per the regress.db ground truth): ~28 reshoot loads across the 3
#      experiments, ~0 many-old→one-cell displacements, every curation count
#      identical pre/post (the dedup drops only the 30 non-survivors' redundant
#      children — all 30 collisions share one image_path, so it is a safe merge).
#   3. Apply for real (back up first — the command prints this reminder):
#        HIMALAYA_DB_PATH=/tmp/upgrade-smoke.db bin/himalaya upgrade-grouping --apply
#   4. Serve the upgraded copy and walk EVERY surface with the browser console
#      open, confirming real data renders with ZERO console errors:
#        HIMALAYA_DB_PATH=/tmp/upgrade-smoke.db bin/himalaya serve --port 8080
#      Surfaces to walk:
#        - Experiment corpus  (loads + samples render under each load)
#        - Experiment grouping-review (fold tree, reshoot flags, structural edits)
#        - Experiment config  (geometry ledger shows derived flight_path / energy
#                              with provenance chips; acquisition timeline renders)
#        - Legacy /samples contact sheet  (thumbnails, phase calls)
#        - Legacy /sample/:id Focus  (trace plate, peaks, assignment rail)
#        - Legacy /series  (folio/scoping/builder overlays)
#        - Legacy /loupe (Inspect)  (detector image, thumbnail gallery, metadata)
#   5. Idempotency by hand: re-run `upgrade-grouping --apply` — the summary must
#      show loads/samples unchanged (loads dedup, owners already claimed). A
#      manual rescan via POST /api/experiments/{id}/scan must be a clean no-op
#      (every derived cell already exists → no new/phantom samples). `serve()`
#      arms no schedulers, so nothing auto-fires while you walk.
# ===========================================================================

using Test
using SQLite, DBInterface, Tables
using Dates
using HimalayaUI

# Shared fixture writers (write_prp / write_setup_info) live in
# test_ingestion_core.jl. Pull them in when running this file standalone.
isdefined(@__MODULE__, :write_prp) || include("test_ingestion_core.jl")

# ---------------------------------------------------------------------------
# count_curation — invariant probe for cross-task tests
# ---------------------------------------------------------------------------

"""
    count_curation(db) -> NamedTuple

COUNT(*) for every table that holds curation / derived state that a migration
or regroup must preserve. The cross-task invariant (Tasks 1/2/4): these counts
are identical pre/post a regroup, and pre/post a dedup minus only the
deduplicated non-survivors' redundant children. Mind the plurals — these are
the EXACT table names.
"""
function count_curation(db::SQLite.DB)
    n(table) = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM $table"))).c
    (auto_peaks            = n("auto_peaks"),
     peak_curations        = n("peak_curations"),
     indices               = n("indices"),
     index_peaks           = n("index_peaks"),
     index_group_members   = n("index_group_members"),
     index_groups          = n("index_groups"),
     assignments           = n("assignments"),
     assignment_members    = n("assignment_members"),
     exposure_tags         = n("exposure_tags"),
     exposure_sources      = n("exposure_sources"),
     series_members        = n("series_members"),
     comparison_members    = n("comparison_members"))
end

# ---------------------------------------------------------------------------
# On-disk fixture writers — give later derivation tasks (Task 2) real files
# ---------------------------------------------------------------------------

"""
    write_stem_fixtures!(data_dir, analysis_dir, stem; horizontal_position_mm, timestamp, kwargs...)

Lay down the complete on-disk fixture set for one exposure stem so a later
scan/regroup can re-derive its load/slot from disk:

  * `<data_dir>/<stem>.tif`  — scan_directory's enumeration key (must exist)
  * `<data_dir>/<stem>.prp`  — carries horizontal_position + timestamp (the
                               partition controls) and the geometry fields
  * `<analysis_dir>/<stem>.dat` — the integration trace input

Control `horizontal_position_mm` (slot partition) and `timestamp` (load
partition: a large inter-frame gap forces a new load) to force a specific
(load, slot) shape. Extra kwargs pass through to `write_prp`.

The caller is responsible for writing ONE shared `setup_info_*.txt` into
`analysis_dir` (via `write_setup_dir!`) — it is per-experiment, not per-stem.
"""
function write_stem_fixtures!(data_dir::AbstractString, analysis_dir::AbstractString,
        stem::AbstractString; horizontal_position_mm::Real,
        timestamp::Union{AbstractString,DateTime} = DateTime(2026, 4, 26, 23, 14, 8),
        kwargs...)
    ts = timestamp isa DateTime ? Dates.format(timestamp, "dd u yyyy HH:MM:SS") : timestamp
    write_prp(joinpath(data_dir, "$stem.prp");
        timestamp = ts, horizontal_position_mm = float(horizontal_position_mm), kwargs...)
    write(joinpath(data_dir, "$stem.tif"), "fake tif")
    write(joinpath(analysis_dir, "$stem.dat"), "fake dat")
    stem
end

"Write the single per-experiment setup_info file into `analysis_dir`."
function write_setup_dir!(analysis_dir::AbstractString; kwargs...)
    path = joinpath(analysis_dir, "setup_info_20260425_181705.txt")
    write_setup_info(path; kwargs...)
    path
end

# ---------------------------------------------------------------------------
# build_legacy_db — the faithful pre-rework fixture
# ---------------------------------------------------------------------------

"""
    build_legacy_db(dir) -> (path, info)

Build a pre-ingestion-rework SQLite DB at `<dir>/legacy.db` in the LEGACY shape
(via `create_schema!`, hand-INSERTed old-shape rows — samples with
name+display_name, exposures with sample_id + bare-stem filename + un-backfilled
NULL experiment_id — then `close`). Does NOT migrate — the caller reopens the path and runs `migrate_schema!` to exercise the
upgrade. Also lays down the on-disk fixture triplets (.tif/.prp/.dat +
setup_info) under `<dir>/data` and `<dir>/analysis` so later derivation tasks
have real files; `experiments.data_dir`/`analysis_dir` point at those dirs.

The fixture embeds the `(experiment_id, filename)` COLLISION that Task 1's dedup
must exercise:

  (a) two exposures with the same bare-stem `filename` but distinct `sample_id`
      (legal under the legacy `(sample_id, filename)` key). One carries a
      `trace_hash` + an `auto` index_groups/auto_peaks curation chain; the other
      is empty. Both point at the SAME `image_path` (the real-data dedup case).
      NOTE: in the legacy shape this collision is legal; the swap to the new
      key is Task 1's job, so `build_legacy_db` deliberately leaves the two rows
      un-deduped and the migration is NOT run here.

PARTITION SHAPE WHEN SCANNED WHOLE (important — verified by probe + the
"build_legacy_db derives to ONE cell when scanned whole" testset below):
**these five single-frame stems derive to ONE load / ONE (load, slot) cell.**
The `reshoot_stems` / `onecell_stems` carry suggestive names and timestamps, but
when `scan_directory` reads the WHOLE dir the gap distribution is unimodal (max
gap < `10×median`, so `_segment_loads` returns a single load) AND the frames are
pure single-frame acquisitions (no multi-frame bursts, so `_cluster_slots` can't
distinguish slot spacing from jitter → one slot). So scanned whole there is no
reshoot split and no many-old→one-cell *separation*; all five frames collapse
into one cell, exercised by the Task 2 "1:1 retrofit + curation preserved" test.

The reshoot (two-loads) and many-old→one-cell shapes are reachable only in
ISOLATION (a fixture that genuinely clusters bursts + a bimodal time gap). Task 2
exercises them via the purpose-built bimodal/co-timed `build_focused_db`; Task 4's
end-to-end round-trip composes ALL of them in one scanned-whole fixture via
`build_legacy_db_full`. `build_legacy_db` is the COLLISION/curation harness, not
a reshoot/many-old harness.

Returns the path plus an `info` NamedTuple of the seeded ids/stems for tests.
The `reshoot`/`one_cell` keys name the seeded sample rows (their stems still
collapse to one cell when scanned whole — see above).
"""
function build_legacy_db(dir::AbstractString)
    data_dir     = joinpath(dir, "data");     mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis");  mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)

    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)   # LEGACY shape: samples(name,display_name); exposures(sample_id,filename), NO experiment_id

    # --- experiment (data_dir/analysis_dir NOT NULL → point at fixture dirs) ---
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'legacy',?,?,?)",
        [dir, data_dir, analysis_dir])

    # ----- (a) collision: two samples, one shared image_path, one bare-stem filename -----
    # 'JC001'/'JC C04' mirror the legacy machine-name / human-label split.
    collide_stem  = "HA_85_422_S2404_0_001"
    collide_image = joinpath(data_dir, "$collide_stem.tif")
    write_stem_fixtures!(data_dir, analysis_dir, collide_stem;
        horizontal_position_mm = 58.9,
        timestamp = DateTime(2026, 4, 26, 23, 14, 8))

    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (10, 1, 'JC001', 'JC C04')")
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (11, 1, 'JC002', 'JC C05')")
    # survivor: carries a trace_hash + a full auto curation chain
    DBInterface.execute(db, """
        INSERT INTO exposures (id, sample_id, filename, image_path, trace_hash)
        VALUES (100, 10, ?, ?, 'deadbeef')""", [collide_stem, collide_image])
    # non-survivor: same (filename, image_path), different sample, no curation
    DBInterface.execute(db, """
        INSERT INTO exposures (id, sample_id, filename, image_path)
        VALUES (101, 11, ?, ?)""", [collide_stem, collide_image])

    # auto curation chain hanging off exposure 100 (must survive a dedup-merge):
    DBInterface.execute(db,
        "INSERT INTO auto_peaks (id, exposure_id, q, intensity) VALUES (1000, 100, 0.123, 5.0)")
    DBInterface.execute(db,
        "INSERT INTO indices (id, exposure_id, phase, basis, kind) VALUES (2000, 100, 'Pn3m', 0.123, 'auto')")
    DBInterface.execute(db,
        "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (3000, 100, 'auto', 1)")
    DBInterface.execute(db,
        "INSERT INTO index_group_members (group_id, index_id) VALUES (3000, 2000)")
    DBInterface.execute(db,
        "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position) VALUES (2000, 1000, 'auto', 1)")

    # ----- (b) reshoot: ONE legacy sample, two stems, same slot, big time gap -----
    reshoot_stems = ["HA_90_500_S2500_0_001", "HA_90_501_S2501_0_001"]
    write_stem_fixtures!(data_dir, analysis_dir, reshoot_stems[1];
        horizontal_position_mm = 70.85, timestamp = DateTime(2026, 4, 26, 23, 30, 0))
    write_stem_fixtures!(data_dir, analysis_dir, reshoot_stems[2];
        horizontal_position_mm = 70.85, timestamp = DateTime(2026, 4, 27, 1, 30, 0)) # +2h gap → new load
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (20, 1, 'HA090', 'HA90 (S01P02)')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (200, 20, ?, ?)",
        [reshoot_stems[1], joinpath(data_dir, "$(reshoot_stems[1]).tif")])
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (201, 20, ?, ?)",
        [reshoot_stems[2], joinpath(data_dir, "$(reshoot_stems[2]).tif")])

    # ----- (c) many-old → one-cell: two legacy samples, single stems each, same load+slot -----
    onecell_stems = ["HA_60_600_S2600_0_001", "HA_60_601_S2601_0_001"]
    write_stem_fixtures!(data_dir, analysis_dir, onecell_stems[1];
        horizontal_position_mm = 63.49, timestamp = DateTime(2026, 4, 26, 23, 40, 0))
    write_stem_fixtures!(data_dir, analysis_dir, onecell_stems[2];
        horizontal_position_mm = 63.50, timestamp = DateTime(2026, 4, 26, 23, 40, 19)) # +19s → same load/slot
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (30, 1, 'HA060a', 'HA60 (S01P03a)')")
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (31, 1, 'HA060b', 'HA60 (S01P03b)')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (300, 30, ?, ?)",
        [onecell_stems[1], joinpath(data_dir, "$(onecell_stems[1]).tif")])
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (301, 31, ?, ?)",
        [onecell_stems[2], joinpath(data_dir, "$(onecell_stems[2]).tif")])

    SQLite.close(db)

    info = (path = path, dir = dir, data_dir = data_dir, analysis_dir = analysis_dir,
            experiment_id = 1,
            collision = (stems = [collide_stem], sample_ids = [10, 11],
                         exposure_ids = [100, 101], image_path = collide_image,
                         survivor_exposure_id = 100),
            reshoot   = (stems = reshoot_stems, sample_id = 20, exposure_ids = [200, 201]),
            one_cell  = (stems = onecell_stems, sample_ids = [30, 31], exposure_ids = [300, 301]))
    return path, info
end

"Open a path, run f(db), always close (mirrors test_ingestion_schema.jl::with_db)."
function with_migration_db(f, path)
    db = SQLite.DB(path)
    try
        return f(db)
    finally
        SQLite.close(db)
    end
end

# ---------------------------------------------------------------------------
# Smoke: the legacy fixture builds and migrate_schema! upgrades it cleanly.
# (Behavior tests — Task 1 dedup / Task 2 regroup — are added by later tasks
#  against these helpers; this file proves only the harness itself.)
# ---------------------------------------------------------------------------

@testset "migration upgrade harness (Task 0)" begin
    @testset "build_legacy_db emits the legacy shape" begin
        path, info = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            # legacy samples shape: name + display_name both present pre-migration
            scols = lowercase.(String.(getproperty.(Tables.rowtable(
                DBInterface.execute(db, "PRAGMA table_info(samples)")), :name)))
            @test "name" in scols
            @test "display_name" in scols
            # legacy exposures shape: experiment_id is present in create_schema!'s
            # final-shape DDL on this branch, but UN-backfilled (NULL) — the
            # faithful pre-rework state the migration's backfill JOIN repairs.
            ehasexp = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id IS NOT NULL"))).c
            @test ehasexp == 0
            # the collision is present and legal under the legacy key
            dupes = Tables.rowtable(DBInterface.execute(db, """
                SELECT filename, COUNT(*) AS n FROM exposures
                 WHERE filename IS NOT NULL GROUP BY filename HAVING n > 1"""))
            @test length(dupes) == 1
            @test dupes[1].filename == info.collision.stems[1]
            # the two colliding exposures share one image_path
            imgs = String.(getproperty.(Tables.rowtable(DBInterface.execute(db,
                "SELECT image_path FROM exposures WHERE id IN (100, 101)")), :image_path))
            @test length(unique(imgs)) == 1
        end
    end

    @testset "on-disk fixture triplets exist for every stem" begin
        path, info = build_legacy_db(mktempdir())
        all_stems = vcat(info.collision.stems, info.reshoot.stems, info.one_cell.stems)
        for stem in all_stems
            @test isfile(joinpath(info.data_dir, "$stem.tif"))   # scan_directory key
            @test isfile(joinpath(info.data_dir, "$stem.prp"))
            @test isfile(joinpath(info.analysis_dir, "$stem.dat"))
        end
        @test !isempty(filter(f -> startswith(f, "setup_info"), readdir(info.analysis_dir)))
        # exposure.filename equals the bare stem (no extension) — the retrofit's
        # `WHERE filename = stem` depends on this.
        with_migration_db(path) do db
            fnames = String.(getproperty.(Tables.rowtable(DBInterface.execute(db,
                "SELECT filename FROM exposures WHERE filename IS NOT NULL")), :filename))
            @test all(!occursin(".", f) for f in fnames)
        end
    end

    @testset "migrate_schema! upgrades the (deduped) legacy fixture cleanly" begin
        # Drop the non-survivor collision row first so the migration's
        # (experiment_id, filename) preflight passes — Task 1 owns the in-migration
        # dedup; here we only prove the harness round-trips through migrate_schema!.
        path, info = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            DBInterface.execute(db, "DELETE FROM exposures WHERE id = ?",
                [info.collision.exposure_ids[2]])
        end

        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end

        with_migration_db(path) do db
            # loads table installed
            tables = lowercase.(String.(getproperty.(Tables.rowtable(DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='table'")), :name)))
            @test "loads" in tables
            # samples.name collapsed; display_name gone
            scols = lowercase.(String.(getproperty.(Tables.rowtable(
                DBInterface.execute(db, "PRAGMA table_info(samples)")), :name)))
            @test "name" in scols && !("display_name" in scols)
            # exposures.experiment_id added + backfilled from the samples JOIN
            ecols = lowercase.(String.(getproperty.(Tables.rowtable(
                DBInterface.execute(db, "PRAGMA table_info(exposures)")), :name)))
            @test "experiment_id" in ecols
            erows = Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id FROM exposures"))
            @test all(r.experiment_id == info.experiment_id for r in erows)
            # the survivor's curation chain rode through the upgrade intact
            c = count_curation(db)
            @test c.auto_peaks == 1
            @test c.indices == 1
            @test c.index_groups == 1
            @test c.index_group_members == 1
            @test c.index_peaks == 1
            # surviving label is the legacy display_name, not the machine name
            sname = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name FROM samples WHERE id = 10"))).name
            @test sname == "JC C04"
        end
    end

    # -----------------------------------------------------------------------
    # Task 1: in-migration delete-only dedup of (experiment_id, filename)
    # collisions that all share one image_path (the real-data case). The
    # build_legacy_db collision fixture (exposures 100 survivor / 101
    # non-survivor, same image_path) drives (a); (c)/(d) inject their own.
    # -----------------------------------------------------------------------
    @testset "Task 1 dedup — same-image_path collision migrates clean" begin
        path, info = build_legacy_db(mktempdir())

        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end

        with_migration_db(path) do db
            # the non-survivor was deleted; the survivor remains.
            ids = sort(Int.(getproperty.(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE filename = ?",
                [info.collision.stems[1]])), :id)))
            @test ids == [info.collision.survivor_exposure_id]

            # the survivor's full curation chain rode through intact. The
            # non-survivor (exposure 101) held no curation, so the survivor's
            # one-of-each chain is exactly what remains.
            c = count_curation(db)
            @test c.auto_peaks == 1 && c.indices == 1 && c.index_groups == 1
            @test c.index_group_members == 1 && c.index_peaks == 1

            # no orphaned FK rows pointing at the deleted non-survivor.
            nsid = info.collision.exposure_ids[2]
            for (tbl, col) in (("auto_peaks", "exposure_id"),
                               ("peak_curations", "exposure_id"),
                               ("indices", "exposure_id"),
                               ("index_groups", "exposure_id"),
                               ("assignments", "exposure_id"),
                               ("assignment_members", "exposure_id"),
                               ("exposure_tags", "exposure_id"))
                n = first(Tables.rowtable(DBInterface.execute(db,
                    "SELECT COUNT(*) AS c FROM $tbl WHERE $col = ?", [nsid]))).c
                @test n == 0
            end

            # the unique index exists and the migration sentinel is recorded.
            idx = Tables.rowtable(DBInterface.execute(db,
                "SELECT name FROM sqlite_master WHERE type='index' AND name='exposures_unique_filename'"))
            @test length(idx) == 1
        end
    end

    @testset "dedup-emptied parent sample is reaped by regroup (no ghost)" begin
        # Collision (a) empties sample 11: its sole exposure (101) is the dedup
        # non-survivor. After regroup it must NOT linger as a ghost (load_id NULL,
        # 0 exposures) — that row is invisible in load rollups but counted by
        # _experiment_stats and listed by /api/samples (both key on experiment_id).
        path, info = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
            # sample 11 is empty post-dedup, pre-regroup (precondition for the ghost).
            @test first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM exposures WHERE sample_id = 11"))).c == 0
            HimalayaUI.regroup_experiment!(db, 1; dry_run = false, analyze = false)
        end
        with_migration_db(path) do db
            # the emptied parent is gone, not lingering.
            @test isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM samples WHERE id = 11")))
            # and no empty, never-grouped sample remains for the experiment.
            @test first(Tables.rowtable(DBInterface.execute(db,
                """SELECT COUNT(*) AS c FROM samples s
                   WHERE s.experiment_id = 1 AND s.load_id IS NULL AND s.merged_into_id IS NULL
                     AND s.id NOT IN (SELECT DISTINCT sample_id FROM exposures
                                      WHERE sample_id IS NOT NULL)"""))).c == 0
        end
    end

    @testset "Task 1 dedup — non-survivor curation is dropped, survivor kept" begin
        # Make the NON-survivor (lower analyzed priority) carry curation so the
        # dedup must delete real children — proves count_curation drops by
        # exactly the non-survivors' redundant rows, not the survivor's.
        path, info = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            # exposure 101 is the non-survivor (NULL trace_hash). Hang an
            # auto_peak + a peak_curation + an exposure_tag off it.
            DBInterface.execute(db,
                "INSERT INTO auto_peaks (id, exposure_id, q, intensity) VALUES (1101, 101, 0.4, 2.0)")
            DBInterface.execute(db,
                "INSERT INTO peak_curations (id, exposure_id, kind, q) VALUES (1102, 101, 'add', 0.5)")
            DBInterface.execute(db,
                "INSERT INTO exposure_tags (id, exposure_id, key, value) VALUES (1103, 101, 'k', 'v')")
        end

        # Count the affected tables directly (count_curation also probes tables
        # the migration adds later, so it can't run pre-migration).
        cnt(db, t) = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM $t"))).c
        with_migration_db(path) do db
            @test cnt(db, "auto_peaks") == 2      # survivor's 1 + non-survivor's 1
            @test cnt(db, "peak_curations") == 1
            @test cnt(db, "exposure_tags") == 1
        end

        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end

        with_migration_db(path) do db
            c = count_curation(db)
            # the non-survivor's redundant children are gone; the survivor's stay.
            @test c.auto_peaks == 1           # only the survivor's auto_peak
            @test c.peak_curations == 0
            @test c.exposure_tags == 0
            @test c.indices == 1              # survivor's, untouched
            @test c.index_peaks == 1
        end
    end

    @testset "Task 1 dedup — re-running migrate_schema! is a no-op" begin
        path, _ = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end
        after_first = with_migration_db(path) do db
            count_curation(db)
        end
        # second run: sentinel already recorded → whole migration skipped.
        with_migration_db(path) do db
            @test HimalayaUI.migrate_schema!(db) === nothing
        end
        after_second = with_migration_db(path) do db
            count_curation(db)
        end
        @test after_first == after_second
    end

    @testset "Task 1 dedup — differing image_path collision still errors" begin
        path, info = build_legacy_db(mktempdir())
        # Point the non-survivor at a DISTINCT image_path → genuinely-distinct
        # files; the dedup must refuse and keep the hard error().
        with_migration_db(path) do db
            DBInterface.execute(db,
                "UPDATE exposures SET image_path = ? WHERE id = ?",
                [joinpath(info.data_dir, "other_image.tif"), info.collision.exposure_ids[2]])
        end
        err = nothing
        @test_throws ErrorException with_migration_db(path) do db
            try
                HimalayaUI.migrate_schema!(db)
            catch e
                err = e
                rethrow()
            end
        end
        @test occursin(info.collision.stems[1], sprint(showerror, err))
    end

    @testset "Task 1 dedup — orphan (NULL sample_id) still trips the orphan error" begin
        path, info = build_legacy_db(mktempdir())
        # Drop the collision non-survivor (so the dedup path isn't what fires)
        # and add an orphan exposure with no derivable experiment_id.
        with_migration_db(path) do db
            DBInterface.execute(db, "DELETE FROM exposures WHERE id = ?",
                [info.collision.exposure_ids[2]])
            DBInterface.execute(db,
                "INSERT INTO exposures (id, sample_id, filename) VALUES (999, NULL, 'orphan')")
        end
        @test_throws ErrorException with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end
    end

    @testset "count_curation returns the twelve named tables" begin
        path, info = build_legacy_db(mktempdir())
        with_migration_db(path) do db
            # drop the collision non-survivor so the migration preflight passes
            DBInterface.execute(db, "DELETE FROM exposures WHERE id = ?",
                [info.collision.exposure_ids[2]])
            HimalayaUI.migrate_schema!(db)
        end
        with_migration_db(path) do db
            c = count_curation(db)
            @test propertynames(c) == (:auto_peaks, :peak_curations, :indices,
                :index_peaks, :index_group_members, :index_groups, :assignments,
                :assignment_members, :exposure_tags, :exposure_sources,
                :series_members, :comparison_members)
            @test all(v -> v isa Integer && v >= 0, values(c))
        end
    end
end

# ---------------------------------------------------------------------------
# Task 2: regroup_experiment! — apply the auto-grouping partition, retrofit rows.
#
# IMPORTANT (real grouping behavior, verified against group_into_samples):
#   * Load segmentation needs a *bimodal* gap distribution to split — a couple of
#     widely-spaced single frames stay ONE load (the unimodal fallback). To force
#     two loads a fixture must cluster several close frames, then a large gap,
#     then several more close frames.
#   * Slot clustering can't distinguish jitter from slot spacing for pure
#     single-frame round-robin acquisitions, so co-timed single frames collapse
#     to ONE slot (the documented Phase-B "single-frame" limitation).
# The Task-0 build_legacy_db fixtures, scanned together in one dir, therefore
# derive to ONE (load, slot) cell — exercised by the collision/curation test
# below. The reshoot (two-loads) and many-old→one-cell shapes need purpose-built
# bimodal/co-timed fixtures, built here via build_focused_db.
# ---------------------------------------------------------------------------

"Build a legacy DB at <dir>/legacy.db and run migrate_schema! (incl. Task 1 dedup)."
function build_migrated_db(dir::AbstractString)
    path, info = build_legacy_db(dir)
    with_migration_db(path) do db
        HimalayaUI.migrate_schema!(db)
    end
    return path, info
end

"""
    build_focused_db(dir; samples, exposures) -> path

Build a focused legacy DB (then migrate it) from explicit `samples` and
`exposures` specs so the derived (load, slot) partition is deterministic.

`samples`   :: Vector of (id, name) — the legacy sample rows (display_name=name).
`exposures` :: Vector of (id, sample_id, stem, hpos, ts::DateTime) — each writes
               its on-disk .tif/.prp/.dat triplet and a legacy exposures row.
"""
function build_focused_db(dir::AbstractString; samples, exposures)
    data_dir = joinpath(dir, "data"); mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis"); mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)
    for e in exposures
        write_stem_fixtures!(data_dir, analysis_dir, e.stem;
            horizontal_position_mm = e.hpos, timestamp = e.ts)
    end
    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'focus',?,?,?)",
        [dir, data_dir, analysis_dir])
    for s in samples
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (?,1,?,?)",
            [s.id, s.name, s.name])
    end
    for e in exposures
        DBInterface.execute(db,
            "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (?,?,?,?)",
            [e.id, e.sample_id, e.stem, joinpath(data_dir, "$(e.stem).tif")])
    end
    SQLite.close(db)
    with_migration_db(path) do db
        HimalayaUI.migrate_schema!(db)
    end
    return path
end

@testset "Task 2 regroup_experiment!" begin
    @testset "1:1 retrofit + curation preserved (collision survivor)" begin
        # build_legacy_db's 5 stems derive to ONE cell; sample 10 (lowest id,
        # carries the surviving auto curation chain) claims it; 11/20/30/31 are
        # displaced. The survivor's row id + human label + curation must ride through.
        path, info = build_migrated_db(mktempdir())
        before = with_migration_db(path) do db; count_curation(db); end
        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, info.experiment_id)
            @test summary.status == :ok

            # Survivor (sample 10) reused in place: id + label kept, load/slot set.
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name, load_id, slot_index, name_source FROM samples WHERE id = 10")))
            @test row.name == "JC C04"
            @test row.load_id !== missing && row.load_id !== nothing
            @test row.slot_index !== missing && row.slot_index !== nothing
            @test row.name_source == "user"

            # The survivor's curation exposure (100) stayed put with a load_id.
            ex = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id, load_id FROM exposures WHERE id = 100")))
            @test ex.sample_id == 10
            @test ex.load_id !== missing && ex.load_id !== nothing
        end
        # count_curation identical before/after — NO exposure deleted/re-inserted.
        after = with_migration_db(path) do db; count_curation(db); end
        @test before == after
    end

    @testset "reshoot → two loads (earliest cell keeps old id+name, later is new auto)" begin
        # Bimodal gap fixture: load A = 3 close frames, +2h gap, load B = 3 close
        # frames, all at one slot position. One legacy sample (20) spans both.
        base = DateTime(2026, 4, 26, 23, 0, 0)
        exps = NamedTuple[]
        for i in 1:3
            push!(exps, (id = 200 + i, sample_id = 20, stem = "RS_A_$(i)_001",
                         hpos = 70.0, ts = base + Dates.Second(2i)))
        end
        for i in 1:3
            push!(exps, (id = 210 + i, sample_id = 20, stem = "RS_B_$(i)_001",
                         hpos = 70.0, ts = base + Dates.Hour(2) + Dates.Second(2i)))
        end
        path = build_focused_db(mktempdir();
            samples = [(id = 20, name = "HA90 label")], exposures = exps)

        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, 1)
            @test summary.loads_created == 2

            # Frames split across two loads / two samples.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, load_id FROM exposures WHERE experiment_id = 1 ORDER BY id"))
            loadids = unique(Int.(getproperty.(rows, :load_id)))
            sampids = unique(Int.(getproperty.(rows, :sample_id)))
            @test length(loadids) == 2
            @test length(sampids) == 2

            # Sample 20 reused for one cell (lowest id), label kept; the other cell
            # is a brand-new 'auto' row.
            @test 20 in sampids
            old = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name, name_source FROM samples WHERE id = 20")))
            @test old.name == "HA90 label"
            @test old.name_source == "user"
            other = first(setdiff(sampids, [20]))
            osrc = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name_source FROM samples WHERE id = ?", [other]))).name_source
            @test osrc == "auto"

            # get_loads_rollup: each load carries its own three exposures.
            roll = HimalayaUI.get_loads_rollup(db, 1)
            @test length(roll) == 2
            for ld in roll
                @test sum(length(sm.exposures) for sm in ld.samples) == 3
            end
        end
    end

    @testset "many-old → one-cell: lowest-id owns, sibling displaced+deleted, metadata carried" begin
        # Two legacy samples, one co-timed single frame each → ONE (load, slot) cell.
        ts0 = DateTime(2026, 4, 26, 23, 40, 0)
        exps = [(id = 300, sample_id = 30, stem = "OC_a_001", hpos = 63.0, ts = ts0),
                (id = 301, sample_id = 31, stem = "OC_b_001", hpos = 63.0,
                 ts = ts0 + Dates.Second(19))]
        path = build_focused_db(mktempdir();
            samples = [(id = 30, name = "HA60a"), (id = 31, name = "HA60b")],
            exposures = exps)
        owner_id, displaced_id = 30, 31

        with_migration_db(path) do db
            # Hang metadata on the soon-to-be-displaced sibling (31): notes + a tag
            # whose key the owner lacks (carries) + a tag whose key the owner holds
            # (must dedup, not throw on sample_tags_unique_key).
            DBInterface.execute(db, "UPDATE samples SET notes = 'displaced note' WHERE id = ?", [displaced_id])
            DBInterface.execute(db, "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'only_displaced', 'v1')", [displaced_id])
            DBInterface.execute(db, "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'shared', 'displaced_val')", [displaced_id])
            DBInterface.execute(db, "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'shared', 'owner_val')", [owner_id])
        end

        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, 1)
            @test summary.samples_displaced == 1

            # Both frames now point at the owner (30) under one load.
            exs = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, load_id FROM exposures WHERE id IN (300, 301) ORDER BY id"))
            @test all(e.sample_id == owner_id for e in exs)
            @test exs[1].load_id == exs[2].load_id

            # Displaced sibling (31) deleted.
            @test isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM samples WHERE id = ?", [displaced_id])))

            # Owner inherited the displaced note via COALESCE (it had none).
            note = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT notes FROM samples WHERE id = ?", [owner_id]))).notes
            @test note == "displaced note"

            # Tag carry: 'only_displaced' migrated; 'shared' kept OWNER's value
            # (dedup-DELETE dropped the displaced dup — no UNIQUE throw).
            tags = Dict(String(r.key) => String(r.value) for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT key, value FROM sample_tags WHERE sample_id = ?", [owner_id])))
            @test tags["only_displaced"] == "v1"
            @test tags["shared"] == "owner_val"
            @test isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM sample_tags WHERE sample_id = ?", [displaced_id])))
        end
    end

    @testset "every file-present exposure gets non-NULL load_id + a sample" begin
        path, info = build_migrated_db(mktempdir())
        with_migration_db(path) do db
            HimalayaUI.regroup_experiment!(db, info.experiment_id)
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, load_id FROM exposures WHERE experiment_id = ?",
                [info.experiment_id]))
            for r in rows
                @test r.load_id !== missing && r.load_id !== nothing
                @test r.sample_id !== missing && r.sample_id !== nothing
            end
        end
    end

    @testset "geometry populated + discrepancies surfaced in the summary" begin
        path, info = build_migrated_db(mktempdir())
        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, info.experiment_id)
            @test haskey(summary, :discrepancies)
            @test summary.discrepancies isa AbstractVector
            @test summary.geometry !== nothing
            g = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT flight_path_m, flight_path_m_source FROM experiments WHERE id = ?",
                [info.experiment_id])))
            @test g.flight_path_m !== missing && g.flight_path_m !== nothing
            @test g.flight_path_m_source != "default"
        end
    end

    @testset "re-running regroup_experiment! is a no-op" begin
        path, info = build_migrated_db(mktempdir())
        snap(db) = (cur = count_curation(db),
            loads = Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [info.experiment_id]))).c),
            samples = Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [info.experiment_id]))).c),
            rollup = HimalayaUI.get_loads_rollup(db, info.experiment_id))
        flatten(roll) = sort([(e.id, sm.sample_id, ld.load_id)
            for ld in roll for sm in ld.samples for e in sm.exposures])

        with_migration_db(path) do db; HimalayaUI.regroup_experiment!(db, info.experiment_id); end
        s1 = with_migration_db(path) do db; snap(db); end
        with_migration_db(path) do db; HimalayaUI.regroup_experiment!(db, info.experiment_id); end
        s2 = with_migration_db(path) do db; snap(db); end

        @test s1.cur == s2.cur
        @test s1.loads == s2.loads
        @test s1.samples == s2.samples
        @test flatten(s1.rollup) == flatten(s2.rollup)
    end

    @testset "empty-metas returns :empty and writes nothing" begin
        dir = mktempdir()
        empty_data = joinpath(dir, "empty_data"); mkpath(empty_data)
        empty_analysis = joinpath(dir, "empty_analysis"); mkpath(empty_analysis)
        path = joinpath(dir, "empty.db")
        db = SQLite.DB(path)
        HimalayaUI.create_schema!(db)
        DBInterface.execute(db,
            "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'empty',?,?,?)",
            [dir, empty_data, empty_analysis])
        SQLite.close(db)
        with_migration_db(path) do db
            HimalayaUI.migrate_schema!(db)
        end
        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, 1)
            @test summary.status == :empty
            @test Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = 1"))).c) == 0
            @test Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = 1"))).c) == 0
        end
    end
end

# ---------------------------------------------------------------------------
# Task 3: upgrade-grouping CLI command — dry-run / apply / idempotency
# ---------------------------------------------------------------------------

"Snapshot all row-level data for loads/samples/exposures in experiment_id."
function snap_experiment(db::SQLite.DB, experiment_id::Int)
    loads = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, load_index, frame_count FROM loads WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name, load_id, slot_index, name_source FROM samples WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))
    exposures = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, sample_id, load_id, filename FROM exposures WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))
    (loads = loads, samples = samples, exposures = exposures,
     curation = count_curation(db))
end

"Run cli_upgrade_grouping with HIMALAYA_DB_PATH=path, capturing stdout to a String."
function run_upgrade(path, args)
    tmpf = tempname()
    open(tmpf, "w") do io
        redirect_stdout(io) do
            withenv("HIMALAYA_DB_PATH" => path) do
                HimalayaUI.cli_upgrade_grouping(args)
            end
        end
    end
    out = read(tmpf, String)
    rm(tmpf; force = true)
    out
end

@testset "Task 3 upgrade-grouping CLI" begin
    @testset "dry-run leaves data unchanged (regroup_experiment! kwarg proof)" begin
        # Unit-level proof: call regroup_experiment! with dry_run=true directly,
        # assert the DB state is row-identical to before.
        path, info = build_migrated_db(mktempdir())
        before = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        with_migration_db(path) do db
            summary = HimalayaUI.regroup_experiment!(db, info.experiment_id; dry_run = true)
            # Summary is populated (what would have happened).
            @test summary.status == :ok
            @test summary.loads_created > 0 || summary.samples_retrofitted > 0 ||
                  summary.samples_created > 0
        end

        after = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        # Every row-level field must be identical — nothing was committed.
        # Full row-state equality (the brief's "byte-identical / full row-state"
        # bar): snap_experiment returns Vector{NamedTuple} for each table. Use
        # `isequal` not `==` — rows hold SQLite NULLs as `missing`, and `==`
        # propagates `missing` (three-valued logic) for any NULL field; `isequal`
        # treats `missing === missing` as true and returns a real Bool. This
        # proves no column of any loads/samples/exposures row drifted under the
        # dry-run rollback.
        @test isequal(before.loads, after.loads)
        @test isequal(before.samples, after.samples)
        @test isequal(before.exposures, after.exposures)
        # No loads should have been created (dry-run rolled back).
        @test isempty(after.loads)
        # Curation untouched.
        @test before.curation == after.curation
    end

    @testset "cli_upgrade_grouping dry-run (default) leaves DB unchanged" begin
        # Exercise arg-parsing + multi-experiment loop via cli_upgrade_grouping
        # called in dry-run mode (no --apply flag).
        path, info = build_migrated_db(mktempdir())
        before = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        out = run_upgrade(path, String[])   # no --apply → dry-run

        after = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        @test isempty(after.loads)      # no loads committed
        @test length(before.samples) == length(after.samples)
        @test length(before.exposures) == length(after.exposures)
        @test before.curation == after.curation
        # Output mentions dry-run.
        @test occursin("DRY-RUN", out) || occursin("dry-run", out)
    end

    @testset "--apply produces the Task 2 end state (loads > 0, partition applied)" begin
        path, info = build_migrated_db(mktempdir())

        out = run_upgrade(path, ["--apply", "--experiment", string(info.experiment_id)])

        with_migration_db(path) do db
            # loads created.
            n_loads = Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?",
                [info.experiment_id]))).c)
            @test n_loads > 0

            # Every exposure has a load_id + sample_id.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, load_id FROM exposures WHERE experiment_id = ?",
                [info.experiment_id]))
            for r in rows
                @test r.load_id !== missing && r.load_id !== nothing
                @test r.sample_id !== missing && r.sample_id !== nothing
            end

            # Output mentions "apply" or status=ok.
            @test occursin("APPLY", out) || occursin("ok", out)
        end
    end

    @testset "--apply twice is a no-op (idempotency)" begin
        path, info = build_migrated_db(mktempdir())
        exp_arg = ["--apply", "--experiment", string(info.experiment_id)]

        # First apply.
        run_upgrade(path, exp_arg)
        snap1 = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        # Second apply — must produce identical state.
        run_upgrade(path, exp_arg)
        snap2 = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        @test length(snap1.loads) == length(snap2.loads)
        @test length(snap1.samples) == length(snap2.samples)
        @test length(snap1.exposures) == length(snap2.exposures)
        @test snap1.curation == snap2.curation

        # Every exposure still has the same sample_id and load_id after the second run.
        s1map = Dict(Int(e.id) => (Int(e.sample_id), Int(e.load_id)) for e in snap1.exposures)
        s2map = Dict(Int(e.id) => (Int(e.sample_id), Int(e.load_id)) for e in snap2.exposures)
        @test s1map == s2map
    end

    @testset "cli_upgrade_grouping --experiment resolves by name" begin
        # Smoke: --experiment <name> is accepted (arg-parse + _resolve_experiment).
        path, info = build_migrated_db(mktempdir())
        out = run_upgrade(path, ["--experiment", "legacy"])
        # Should not error; dry-run output produced.
        @test occursin("DRY-RUN", out) || occursin("dry-run", out)
    end

    @testset "unreachable data_dir skips with a message" begin
        path, info = build_migrated_db(mktempdir())
        # Point the experiment at a non-existent dir.
        with_migration_db(path) do db
            DBInterface.execute(db,
                "UPDATE experiments SET data_dir = ? WHERE id = ?",
                ["/nonexistent/path/$(rand(1:999999))", info.experiment_id])
        end

        out = run_upgrade(path, ["--apply", "--experiment", string(info.experiment_id)])

        # No loads created — the skip path fired.
        with_migration_db(path) do db
            n = Int(first(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?",
                [info.experiment_id]))).c)
            @test n == 0
        end
        @test occursin("SKIP", out)
    end
end

# ---------------------------------------------------------------------------
# Task 4: end-to-end round-trip — collisions + curation + reshoot + many-old →
# one-cell, ALL composed in one fixture that GENUINELY produces those shapes
# when the CLI scans its data_dir whole. (build_legacy_db can't: scanned whole
# its 5 single-frame stems collapse to one cell — see the carried-finding
# discussion in build_legacy_db's docstring. So Task 4 owns build_legacy_db_full,
# whose bimodal time gap + multi-frame bursts make the derivation produce two
# loads with a reshoot split and a many-old→one-cell collapse.)
#
# Asserted partition (verified before any migration assertion below):
#   Load 1, slot 1 (pos 60.0): OC_a + OC_b + collision-survivor stem  (3 frames)
#   Load 1, slot 2 (pos 70.0): reshoot batch 1                        (2 frames)
#   Load 2, slot 1 (pos 70.0): reshoot batch 2                        (2 frames)
# → 2 loads, 3 cells.  Sample assignment in build_legacy_db_full:
#   - collision: sample 10 (survivor, curation chain) + 11 (non-survivor, same
#     image_path) share filename COLLIDE; the survivor's frame sits in cell(1,1).
#   - many-old → one-cell: samples 30/31 own OC_a/OC_b → both land in cell(1,1).
#     With the collision survivor that is THREE old samples in one cell; lowest
#     id (10) owns, 30 & 31 are displaced+deleted (metadata carried).
#   - reshoot: sample 20 owns BOTH reshoot batches → cell(1,2) keeps id 20+name,
#     cell(2,1) becomes a new 'auto' sample.
# ---------------------------------------------------------------------------

"""
    build_legacy_db_full(dir) -> (path, info)

Like `build_legacy_db`, but the on-disk fixture is engineered so that when
`scan_directory` reads the WHOLE data_dir the derivation produces the full set
of migration shapes at once (bimodal time gap → two loads; multi-frame bursts →
≥2 slots; a reshoot specimen split across loads; three old samples collapsing
into one cell). See the block comment above for the exact asserted partition.

The collision (same `image_path`, two distinct `sample_id`) is folded INTO the
many-old cell: the surviving collision exposure (id 100, carrying the curation
chain) is one of cell(1,1)'s three frames, so after Task 1's dedup + Task 4's
regroup, sample 10 owns cell(1,1) and carries both the survived curation and the
absorbed siblings' metadata.

`info` mirrors build_legacy_db's shape plus a `nonsurvivor_curation` count of the
curation children hung off the collision non-survivor (exposure 101) — the exact
amount `count_curation` must drop across the dedup.
"""
function build_legacy_db_full(dir::AbstractString)
    data_dir     = joinpath(dir, "data");      mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis");   mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)

    base = DateTime(2026, 4, 26, 23, 0, 0)
    S(n) = Dates.Second(n)

    collide_stem  = "HA_85_422_S2404_0_001"
    collide_image = joinpath(data_dir, "$collide_stem.tif")

    # (stem, hpos, timestamp). One on-disk frame each (the collision is two
    # exposure ROWS sharing collide_stem's single on-disk file).
    onecell_stems = ["OC_a_S2600_0_001", "OC_b_S2601_0_001"]
    reshoot_stems = ["RS_S2500_0_001", "RS_S2500_0_002", "RS_S2501_0_001", "RS_S2501_0_002"]
    disk = [
        # Load 1, slot 1 (pos 60.0): many-old trio (OC_a, OC_b, collision survivor)
        (onecell_stems[1], 60.0, base + S(0)),
        (onecell_stems[2], 60.0, base + S(20)),
        (collide_stem,     60.0, base + S(25)),
        # Load 1, slot 2 (pos 70.0): reshoot batch 1 (burst)
        (reshoot_stems[1], 70.0, base + S(40)),
        (reshoot_stems[2], 70.0, base + S(42)),
        # Load 2 (pos 70.0): reshoot batch 2 (burst), +2h gap
        (reshoot_stems[3], 70.0, base + Dates.Hour(2) + S(0)),
        (reshoot_stems[4], 70.0, base + Dates.Hour(2) + S(2)),
    ]
    for (stem, hpos, ts) in disk
        write_stem_fixtures!(data_dir, analysis_dir, stem;
            horizontal_position_mm = hpos, timestamp = ts)
    end

    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'legacyfull',?,?,?)",
        [dir, data_dir, analysis_dir])

    # samples: 10 collision survivor (human label), 11 collision non-survivor,
    # 20 reshoot specimen, 30/31 many-old pair.
    for (id, name, dname) in [(10,"JC001","JC C04"), (11,"JC002","JC C05"),
                              (20,"HA090","HA90 (S01P02)"),
                              (30,"HA060a","HA60a"), (31,"HA060b","HA60b")]
        DBInterface.execute(db,
            "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (?,1,?,?)",
            [id, name, dname])
    end

    # collision survivor (100, sample 10) + a full auto curation chain.
    DBInterface.execute(db, """
        INSERT INTO exposures (id, sample_id, filename, image_path, trace_hash)
        VALUES (100, 10, ?, ?, 'deadbeef')""", [collide_stem, collide_image])
    DBInterface.execute(db,
        "INSERT INTO auto_peaks (id, exposure_id, q, intensity) VALUES (1000, 100, 0.123, 5.0)")
    DBInterface.execute(db,
        "INSERT INTO indices (id, exposure_id, phase, basis, kind) VALUES (2000, 100, 'Pn3m', 0.123, 'auto')")
    DBInterface.execute(db,
        "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (3000, 100, 'auto', 1)")
    DBInterface.execute(db,
        "INSERT INTO index_group_members (group_id, index_id) VALUES (3000, 2000)")
    DBInterface.execute(db,
        "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position) VALUES (2000, 1000, 'auto', 1)")

    # collision NON-survivor (101, sample 11): same (filename, image_path). Hang
    # curation on it so the dedup must drop REAL children (a meaningful
    # count_curation delta, not a vacuous pre==post). nonsurvivor_curation tracks
    # the exact drop: 1 auto_peak + 1 peak_curation + 1 exposure_tag = 3 rows.
    DBInterface.execute(db, """
        INSERT INTO exposures (id, sample_id, filename, image_path)
        VALUES (101, 11, ?, ?)""", [collide_stem, collide_image])
    DBInterface.execute(db,
        "INSERT INTO auto_peaks (id, exposure_id, q, intensity) VALUES (1101, 101, 0.4, 2.0)")
    DBInterface.execute(db,
        "INSERT INTO peak_curations (id, exposure_id, kind, q) VALUES (1102, 101, 'add', 0.5)")
    DBInterface.execute(db,
        "INSERT INTO exposure_tags (id, exposure_id, key, value) VALUES (1103, 101, 'k', 'v')")
    nonsurvivor_curation = (auto_peaks = 1, peak_curations = 1, exposure_tags = 1)

    # reshoot specimen 20: all four reshoot frames.
    for (i, stem) in enumerate(reshoot_stems)
        DBInterface.execute(db,
            "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (?, 20, ?, ?)",
            [200 + i, stem, joinpath(data_dir, "$stem.tif")])
    end
    # many-old pair: 30 owns OC_a, 31 owns OC_b.
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (300, 30, ?, ?)",
        [onecell_stems[1], joinpath(data_dir, "$(onecell_stems[1]).tif")])
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (301, 31, ?, ?)",
        [onecell_stems[2], joinpath(data_dir, "$(onecell_stems[2]).tif")])

    SQLite.close(db)

    info = (path = path, dir = dir, data_dir = data_dir, analysis_dir = analysis_dir,
            experiment_id = 1,
            collide_stem = collide_stem,
            survivor_exposure_id = 100, nonsurvivor_exposure_id = 101,
            survivor_sample_id = 10, nonsurvivor_sample_id = 11,
            reshoot_sample_id = 20, reshoot_stems = reshoot_stems,
            onecell_sample_ids = [30, 31], onecell_stems = onecell_stems,
            nonsurvivor_curation = nonsurvivor_curation)
    return path, info
end


@testset "Task 4 end-to-end round-trip" begin
    @testset "fixture derives the asserted multi-shape partition when scanned whole" begin
        # FIRST prove the fixture shape (not assume it): run the pure derivation
        # the CLI runs, and assert two loads / three cells with the exact frame
        # membership. Everything downstream depends on this partition.
        path, info = build_legacy_db_full(mktempdir())
        metas = HimalayaUI.scan_directory(info.data_dir, info.analysis_dir)
        @test length(metas) == 7
        result = HimalayaUI.group_into_samples(metas)
        @test length(result.loads) == 2

        l1, l2 = result.loads[1], result.loads[2]
        @test l1.frame_count == 5
        @test l2.frame_count == 2
        @test length(l1.samples) == 2          # slot 1 (many-old trio) + slot 2 (reshoot b1)
        @test length(l2.samples) == 1          # reshoot b2

        # cell(1,1) holds the many-old trio (OC_a, OC_b, collision survivor stem)
        cell11 = Set(e.stem for e in l1.samples[1].exposures)
        @test cell11 == Set([info.onecell_stems[1], info.onecell_stems[2], info.collide_stem])
        # cell(1,2) + cell(2,1) hold the reshoot's two bursts.
        cell12 = Set(e.stem for e in l1.samples[2].exposures)
        cell21 = Set(e.stem for e in l2.samples[1].exposures)
        @test cell12 == Set(info.reshoot_stems[1:2])
        @test cell21 == Set(info.reshoot_stems[3:4])
    end

    @testset "build → migrate_schema! (schema + dedup) → upgrade-grouping --apply" begin
        path, info = build_legacy_db_full(mktempdir())

        # curation BEFORE the migration (probe the affected tables directly —
        # count_curation also reads tables migrate_schema! adds, so it can't run
        # pre-migration). Both the survivor (100) and non-survivor (101) carry
        # a one-each chain at this point.
        cnt(db, t) = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM $t"))).c)
        before_raw = with_migration_db(path) do db
            (auto_peaks = cnt(db, "auto_peaks"),
             peak_curations = cnt(db, "peak_curations"),
             exposure_tags = cnt(db, "exposure_tags"),
             indices = cnt(db, "indices"),
             index_peaks = cnt(db, "index_peaks"))
        end
        @test before_raw.auto_peaks == 2          # survivor 1 + non-survivor 1
        @test before_raw.peak_curations == 1      # non-survivor's
        @test before_raw.exposure_tags == 1       # non-survivor's

        # --- schema half: migrate_schema! upgrades + runs Task 1's dedup ---
        with_migration_db(path) do db
            @test HimalayaUI.migrate_schema!(db) === nothing  # opens/migrates clean
        end

        # curation AFTER the dedup but BEFORE regroup.
        after_dedup = with_migration_db(path) do db; count_curation(db); end
        # The dedup dropped exactly the non-survivor's redundant children.
        @test after_dedup.auto_peaks     == before_raw.auto_peaks - info.nonsurvivor_curation.auto_peaks
        @test after_dedup.peak_curations == before_raw.peak_curations - info.nonsurvivor_curation.peak_curations
        @test after_dedup.exposure_tags  == before_raw.exposure_tags - info.nonsurvivor_curation.exposure_tags
        # The survivor's chain is intact.
        @test after_dedup.auto_peaks == 1 && after_dedup.indices == 1 && after_dedup.index_peaks == 1
        # The collision is deduped to the survivor only.
        with_migration_db(path) do db
            ids = sort(Int.(getproperty.(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM exposures WHERE filename = ?", [info.collide_stem])), :id)))
            @test ids == [info.survivor_exposure_id]
        end

        # --- data half: exercise the CLI --apply path (Task 3 wrapper) ---
        out = run_upgrade(path, ["--apply", "--experiment", string(info.experiment_id)])
        @test occursin("APPLY", out)
        @test occursin("ok", out)

        with_migration_db(path) do db
            # loads > 0 with correct frame_counts (5 / 2 from the partition).
            loads = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, load_index, frame_count FROM loads WHERE experiment_id = ? ORDER BY load_index",
                [info.experiment_id]))
            @test length(loads) == 2
            @test loads[1].frame_count == 5
            @test loads[2].frame_count == 2

            # Partition: one sample per derived cell → 3 live cells under 2 loads.
            roll = HimalayaUI.get_loads_rollup(db, info.experiment_id)
            @test length(roll) == 2
            cells_per_load = sort([length(ld.samples) for ld in roll])
            @test cells_per_load == [1, 2]

            # Reshoot: sample 20 keeps the EARLIEST cell (id + human name); the
            # later cell is a NEW 'auto' sample.
            r20 = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name, name_source, load_id, slot_index FROM samples WHERE id = 20")))
            @test r20.name == "HA90 (S01P02)"
            @test r20.name_source == "user"
            # The reshoot's later cell is a fresh auto sample (id ∉ the seeded set).
            seeded = Set([10, 11, 20, 30, 31])
            reshoot_stem_sids = Set(Int(r.sample_id) for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT DISTINCT sample_id FROM exposures WHERE filename IN (?, ?)",
                [info.reshoot_stems[3], info.reshoot_stems[4]])))
            @test 20 ∉ reshoot_stem_sids                        # later batch moved off 20
            new_reshoot = first(setdiff(reshoot_stem_sids, seeded))
            nr = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name_source FROM samples WHERE id = ?", [new_reshoot])))
            @test nr.name_source == "auto"

            # Many-old → one cell: lowest id (10) owns cell(1,1); 30 & 31 deleted.
            @test isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM samples WHERE id IN (30, 31)")))
            owner_sids = Set(Int(r.sample_id) for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT DISTINCT sample_id FROM exposures WHERE filename IN (?, ?, ?)",
                [info.onecell_stems[1], info.onecell_stems[2], info.collide_stem])))
            @test owner_sids == Set([info.survivor_sample_id])  # all three frames now on sample 10
            @test first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name FROM samples WHERE id = 10"))).name == "JC C04"  # human label kept

            # Every file-present exposure relinked: non-NULL load_id + a sample.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, sample_id, load_id FROM exposures WHERE experiment_id = ?",
                [info.experiment_id]))
            for r in rows
                @test r.load_id   !== missing && r.load_id   !== nothing
                @test r.sample_id !== missing && r.sample_id !== nothing
            end

            # count_curation post-regroup == post-dedup (regroup deletes/re-inserts
            # NO exposure → exposure-keyed curation is untouched by the data half).
            @test count_curation(db) == after_dedup

            # Geometry derived from PRP/setup (flight_path from setup_info).
            g = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT flight_path_m, flight_path_m_source, energy_kev, energy_kev_source, " *
                "q_units, q_units_source FROM experiments WHERE id = ?", [info.experiment_id])))
            @test g.flight_path_m !== missing && g.flight_path_m !== nothing
            @test g.flight_path_m_source != "default"
            @test g.energy_kev !== missing && g.energy_kev !== nothing
            @test g.energy_kev_source != "default"
            # q_units has no writer → stays NULL/'default' (DO NOT assert populated).
            @test g.q_units === missing || g.q_units === nothing
            @test g.q_units_source == "default"
        end
    end

    @testset "summary reshoot/displaced counts match the fixture's known shape" begin
        # Record the summary counts so a future real-dev-db run can be eyeballed
        # against ≈28 reshoots / ≈0 many-old→one-cell. On THIS fixture: 2 loads
        # (one reshoot split), 3 cells, 2 retrofitted (10 + 20), 1 created (the
        # reshoot later cell), 2 displaced (30 + 31 → the many-old collapse).
        # `reshoots` = old samples whose stems span ≥2 cells; here ONLY sample 20
        # (its two bursts fall in cell(1,2) + cell(2,1)) → 1.
        path, info = build_legacy_db_full(mktempdir())
        with_migration_db(path) do db; HimalayaUI.migrate_schema!(db); end
        summary = with_migration_db(path) do db
            HimalayaUI.regroup_experiment!(db, info.experiment_id)
        end
        @test summary.status == :ok
        @test summary.loads_created == 2
        @test summary.reshoots == 1            # sample 20 spans 2 cells (the only split)
        @test summary.cells == 3
        @test summary.samples_retrofitted == 2
        @test summary.samples_created == 1     # the reshoot's later cell
        @test summary.samples_displaced == 2   # many-old collapse (30, 31)
        @test summary.exposures_no_file == 0
        # Every on-disk frame relinked (7 exposures: 1 collision survivor + 4
        # reshoot + 2 many-old; the non-survivor was deduped away pre-regroup).
        @test summary.exposures_relinked == 7
    end

    @testset "a second --apply is a clean no-op (idempotency)" begin
        path, info = build_legacy_db_full(mktempdir())
        with_migration_db(path) do db; HimalayaUI.migrate_schema!(db); end
        exp_arg = ["--apply", "--experiment", string(info.experiment_id)]

        run_upgrade(path, exp_arg)
        snap1 = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        run_upgrade(path, exp_arg)
        snap2 = with_migration_db(path) do db; snap_experiment(db, info.experiment_id); end

        # Full row-state identity across loads/exposures + curation.
        # `isequal` (not `==`) for the row vectors — NULL columns surface as
        # `missing` and `==` would propagate `missing` instead of a Bool.
        @test isequal(snap1.loads, snap2.loads)
        @test isequal(snap1.exposures, snap2.exposures)
        @test snap1.curation == snap2.curation

        # Samples: full row identity, including name_source. The capstone originally
        # caught a name_source idempotency defect here (a second --apply flipped a
        # regroup-auto-created reshoot cell 'auto'→'user' because the retrofit branch
        # stamped 'user' unconditionally). Fixed in regroup_experiment! (ingest.jl):
        # the stamp is now `CASE WHEN load_id IS NULL THEN 'user' ELSE name_source END`
        # — 'user' only on a genuine first adoption of a pre-rework manifest sample
        # (load_id NULL); re-claiming a sample regroup already placed leaves it 'auto'.
        @test isequal(snap1.samples, snap2.samples)

        # cheap_change_check confirms no-op (false): every on-disk image is already
        # persisted. NOTE: post-dedup `persisted` can differ from the on-disk .tif
        # count only if data_dir holds stray un-persisted .tif (e.g. calibration
        # frames); this fixture writes none, so on_disk == persisted → false.
        with_migration_db(path) do db
            @test HimalayaUI.cheap_change_check(db, info.experiment_id) == false
        end
    end
end

# ---------------------------------------------------------------------------
# Real-data shape: manifest-era filename is the on-disk stem with the trailing
# `_<rep>_<frame>` suffix STRIPPED (live dev-db: every curated `exposures.filename`
# = scan stem minus `_0_001`). scan_directory keeps the suffix, so a relink that
# matches `filename == scan_stem` finds NOTHING and orphans all curation.
# build_legacy_db's fixtures hide this by setting filename = the full disk stem;
# this fixture sets filename to the TRUNCATED form, the way production stores it.
# ---------------------------------------------------------------------------
@testset "regroup relinks manifest-era truncated filenames (P0-a)" begin
    dir          = mktempdir()
    data_dir     = joinpath(dir, "data");      mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis");  mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)

    # On disk: the FULL stem (what scan_directory enumerates).
    # In the DB:  the TRUNCATED stem (what the manifest pipeline stored).
    full_stem   = "JC_C01_1_S2453_0_001"
    db_filename = "JC_C01_1_S2453"
    write_stem_fixtures!(data_dir, analysis_dir, full_stem;
        horizontal_position_mm = 58.9, timestamp = DateTime(2026, 4, 26, 23, 14, 8))

    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'legacy',?,?,?)",
        [dir, data_dir, analysis_dir])
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (10, 1, 'JC001', 'JC C01')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path, trace_hash) VALUES (100, 10, ?, ?, 'deadbeef')",
        [db_filename, joinpath(data_dir, "$full_stem.tif")])
    # Curation hanging off the exposure — must survive the relink.
    DBInterface.execute(db,
        "INSERT INTO auto_peaks (id, exposure_id, q, intensity) VALUES (1000, 100, 0.123, 5.0)")
    SQLite.close(db)

    with_migration_db(path) do db
        HimalayaUI.migrate_schema!(db)
    end

    summary = with_migration_db(path) do db
        HimalayaUI.regroup_experiment!(db, 1; dry_run = false)
    end

    # The curated exposure must be RELINKED (retrofit in place), not orphaned.
    @test summary.exposures_relinked == 1
    @test summary.samples_retrofitted == 1

    with_migration_db(path) do db
        # Same exposure id, with its curation still attached.
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, sample_id, load_id FROM exposures WHERE id = 100")))
        @test row.sample_id == 10           # retrofitted onto the existing sample
        @test row.load_id !== missing       # placed under a load
        # filename rewritten to the full scan stem so future rescans dedup correctly.
        @test row.filename == full_stem
        # curation preserved.
        @test first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM auto_peaks WHERE exposure_id = 100"))).c == 1
    end
end

# ---------------------------------------------------------------------------
# Ingest-everything: a scanned file with NO pre-rework DB row (an un-manifested
# sample — the live dev-db has hundreds: HA_/RY_ series + calibration) must be
# INSERTED as a real exposure so it surfaces as a sample, not minted as an empty
# sample row. (Decision 2026-06-20: experiment == directory, ingest everything.)
# ---------------------------------------------------------------------------
@testset "regroup inserts exposures for un-manifested scanned files (ingest-everything)" begin
    dir          = mktempdir()
    data_dir     = joinpath(dir, "data");      mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis");  mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)

    # One curated file (DB row, truncated filename) + one un-manifested file (on
    # disk only). Distinct horizontal positions → two slots, one load.
    curated_full = "JC_C01_1_S2453_0_001"
    curated_db   = "JC_C01_1_S2453"
    newfile_full = "HA_5_010_S1965_0_001"   # no DB row — a fresh-scan-only file
    write_stem_fixtures!(data_dir, analysis_dir, curated_full;
        horizontal_position_mm = 58.9, timestamp = DateTime(2026, 4, 26, 23, 14, 8))
    write_stem_fixtures!(data_dir, analysis_dir, newfile_full;
        horizontal_position_mm = 70.85, timestamp = DateTime(2026, 4, 26, 23, 14, 20))

    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'legacy',?,?,?)",
        [dir, data_dir, analysis_dir])
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name, display_name) VALUES (10, 1, 'JC001', 'JC C01')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, image_path) VALUES (100, 10, ?, ?)",
        [curated_db, joinpath(data_dir, "$curated_full.tif")])
    SQLite.close(db)

    with_migration_db(path) do db
        HimalayaUI.migrate_schema!(db)
    end

    # analyze=false: fixtures carry fake tif/dat; isolate the insert from peak-finding.
    summary = with_migration_db(path) do db
        HimalayaUI.regroup_experiment!(db, 1; dry_run = false, analyze = false)
    end

    @test summary.exposures_relinked == 1     # curated file relinked
    @test summary.exposures_inserted == 1     # un-manifested file inserted

    with_migration_db(path) do db
        # The new file now has a real exposure row with a usable image_path.
        new = Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, sample_id, load_id, image_path FROM exposures WHERE filename = ?",
            [newfile_full]))
        @test length(new) == 1
        @test new[1].sample_id !== missing
        @test new[1].load_id !== missing
        @test new[1].image_path !== missing   # image route filters WHERE image_path IS NOT NULL
        # No empty samples: every sample owns ≥1 exposure.
        empties = first(Tables.rowtable(DBInterface.execute(db,
            """SELECT COUNT(*) AS c FROM samples s
               WHERE s.experiment_id = 1
                 AND NOT EXISTS (SELECT 1 FROM exposures e WHERE e.sample_id = s.id)"""))).c
        @test empties == 0
    end
end

# ---------------------------------------------------------------------------
# Indistinguishable-from-fresh-scan (Jonathan's bar, 2026-06-20): a migrated
# experiment must look exactly like a freshly-ingested one. Direct proof — run
# scan_and_group! on a fresh DB and regroup_experiment! on a legacy DB built from
# the SAME directory, then compare the persisted loads/samples/exposures.
#
# The bar applies to STRUCTURE (load partition + start/end/frame_count, slot/load
# placement, grouping_source, and every exposure field) — NOT to preserved
# CURATION. A migrated sample keeps its human label (the real dev-db's
# `display_name`, e.g. "2-2 + LL37 1:1", collapses into `name` with
# name_source='user'); a fresh scan auto-derives "JC_C01_1 (S01P01)"/'auto'. That
# difference is the whole point of preserving curation, so name/name_source are
# excluded from the structural comparison and asserted separately below.
# ---------------------------------------------------------------------------
@testset "migrated experiment == fresh scan (structural equivalence)" begin
    # Shared on-disk fixture: a bimodal-gap set → 2 loads, several slots.
    src          = mktempdir()
    data_dir     = joinpath(src, "data");      mkpath(data_dir)
    analysis_dir = joinpath(src, "analysis");  mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)
    stems = ["JC_C01_1_S2453_0_001", "JC_C01_2_S2454_0_001", "JC_C01_3_S2455_0_001",
             "JC_C02_1_S2460_0_001", "JC_C02_2_S2461_0_001"]
    # First three co-timed (one load), +2h gap, last two co-timed (second load).
    base = DateTime(2026, 4, 26, 23, 0, 0)
    ts   = [base, base + Second(20), base + Second(40),
            base + Hour(2), base + Hour(2) + Second(20)]
    hpos = [58.9, 63.5, 70.8, 58.9, 63.5]
    for (i, s) in enumerate(stems)
        write_stem_fixtures!(data_dir, analysis_dir, s;
            horizontal_position_mm = hpos[i], timestamp = ts[i])
    end

    # Snapshot the structural (non-id, non-curation) shape of an experiment.
    function shape(db, eid)
        loads = Tables.rowtable(DBInterface.execute(db,
            "SELECT load_index, frame_count, start_time, end_time FROM loads WHERE experiment_id = ? ORDER BY load_index", [eid]))
        samples = Tables.rowtable(DBInterface.execute(db,
            """SELECT l.load_index, s.slot_index, s.grouping_source
               FROM samples s JOIN loads l ON s.load_id = l.id
               WHERE s.experiment_id = ? ORDER BY l.load_index, s.slot_index""", [eid]))
        exposures = Tables.rowtable(DBInterface.execute(db,
            """SELECT e.filename, l.load_index, s.slot_index, e.timestamp,
                      e.horizontal_position, e.scan_id, e.frame_no
               FROM exposures e JOIN samples s ON e.sample_id = s.id JOIN loads l ON e.load_id = l.id
               WHERE e.experiment_id = ? ORDER BY e.filename""", [eid]))
        (loads = loads, samples = samples, exposures = exposures)
    end

    # (a) Fresh scan.
    fresh_path = joinpath(src, "fresh.db")
    fdb = SQLite.DB(fresh_path); HimalayaUI.create_schema!(fdb)
    DBInterface.execute(fdb,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'fresh',?,?,?)",
        [src, data_dir, analysis_dir])
    SQLite.close(fdb)
    fresh = with_migration_db(fresh_path) do db
        HimalayaUI.migrate_schema!(db)
        HimalayaUI.scan_and_group!(db, 1; analyze = false)
        shape(db, 1)
    end

    # (b) Legacy DB from the same dir: every stem a manifest-era row (truncated
    # filename, one legacy sample), then migrate + regroup.
    leg_path = joinpath(src, "legacy.db")
    ldb = SQLite.DB(leg_path); HimalayaUI.create_schema!(ldb)
    DBInterface.execute(ldb,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'legacy',?,?,?)",
        [src, data_dir, analysis_dir])
    for (i, s) in enumerate(stems)
        sid = 100 + i
        # Distinct grouping_source → the retrofit must re-derive it to match a fresh
        # scan ('auto_position'), not keep the stale legacy value.
        DBInterface.execute(ldb,
            "INSERT INTO samples (id, experiment_id, name, display_name, grouping_source) VALUES (?, 1, ?, ?, 'legacy_manual')",
            [sid, "M$i", "Manifest $i"])
        DBInterface.execute(ldb,
            "INSERT INTO exposures (experiment_id, sample_id, filename, image_path) VALUES (1, ?, ?, ?)",
            [sid, replace(s, r"_\d+_\d+$" => ""), joinpath(data_dir, "$s.tif")])
    end
    SQLite.close(ldb)
    migrated = with_migration_db(leg_path) do db
        HimalayaUI.migrate_schema!(db)
        HimalayaUI.regroup_experiment!(db, 1; dry_run = false, analyze = false)
        shape(db, 1)
    end

    @test isequal(fresh.loads, migrated.loads)
    @test isequal(fresh.samples, migrated.samples)        # structure only (no name)
    @test isequal(fresh.exposures, migrated.exposures)

    # Curation IS preserved (the deliberate, correct divergence): the migrated
    # samples keep their human manifest labels with name_source='user', while the
    # fresh scan auto-derives names with name_source='auto'.
    mig_names = with_migration_db(leg_path) do db
        Tables.rowtable(DBInterface.execute(db,
            """SELECT s.name, s.name_source FROM samples s JOIN loads l ON s.load_id = l.id
               WHERE s.experiment_id = 1 ORDER BY l.load_index, s.slot_index"""))
    end
    @test all(r -> r.name_source == "user", mig_names)
    @test startswith(mig_names[1].name, "Manifest")       # human label kept, not auto-derived
end

# ---------------------------------------------------------------------------
# Production trace naming: a real pre-rework DB's integration traces are the
# per-acquisition `_tot.dat` totals (HA_5_010_S1965_tot.dat). Single-frame ingest
# (the migration resolves the funnel's detected `{name}_0_001.*` pattern from
# disk, dropping the frame suffix) names the exposure by the acquisition stem
# (HA_5_010_S1965) — identical to a fresh scan. The migration's analyze step must
# resolve the shared `_tot.dat` so inserted exposures get peaks, not just exist.
# ---------------------------------------------------------------------------
@testset "migration analyzes inserted exposures via per-acquisition _tot.dat" begin
    dir          = mktempdir()
    data_dir     = joinpath(dir, "data");      mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis");  mkpath(analysis_dir)
    write_setup_dir!(analysis_dir)

    full = "HA_5_010_S1965_0_001"   # on-disk per-frame file
    acq  = "HA_5_010_S1965"         # acquisition stem == ingested exposures.filename
    write_stem_fixtures!(data_dir, analysis_dir, full;
        horizontal_position_mm = 58.9, timestamp = DateTime(2026, 4, 26, 23, 14, 8))
    # Real trace bytes named by the ACQUISITION stem (frame suffix dropped).
    fixture = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(fixture, joinpath(analysis_dir, "HA_5_010_S1965_tot.dat"); force = true)

    path = joinpath(dir, "legacy.db")
    db = SQLite.DB(path)
    HimalayaUI.create_schema!(db)
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir, config) VALUES (1,'legacy',?,?,?,?)",
        [dir, data_dir, analysis_dir, "[files]\nintegration = \"{name}_tot.dat\"\n"])
    SQLite.close(db)

    with_migration_db(path) do db
        HimalayaUI.migrate_schema!(db)
    end
    with_migration_db(path) do db
        HimalayaUI.regroup_experiment!(db, 1; dry_run = false, analyze = true)
    end

    with_migration_db(path) do db
        eid = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE filename = ?", [acq]))).id
        npeaks = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM auto_peaks WHERE exposure_id = ?", [eid]))).c
        @test npeaks > 0   # inserted exposure got analyzed via the frame-suffix fallback
    end
end

# ---------------------------------------------------------------------------
# backfill_load_sessions! — idempotent backfill of loads.session_id
# ---------------------------------------------------------------------------
@testset "backfill_load_sessions! populates session_id from start_times" begin
    # Build a fresh open_db so migrate_schema! has already run (incl. the backfill
    # sentinel). Then directly test backfill_load_sessions! on a DB that has loads
    # with NULL session_id.
    path = joinpath(mktempdir(), "backfill.db")
    db = HimalayaUI.open_db(path)

    # Insert one experiment with 3 loads: two within 1 h (session 1), then
    # a >3 h gap, then one more (session 2).
    DBInterface.execute(db,
        "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1,'bftest','/t','/d','/a')")
    t0 = "2026-04-25T10:00:00"
    t1 = "2026-04-25T10:30:00"
    t2 = "2026-04-25T14:00:00"   # +4 h gap → new session
    # Insert loads with NULL session_id directly (bypassing create_load! which sets it).
    DBInterface.execute(db,
        "INSERT INTO loads (id, experiment_id, load_index, start_time) VALUES (101, 1, 1, ?)", [t0])
    DBInterface.execute(db,
        "INSERT INTO loads (id, experiment_id, load_index, start_time) VALUES (102, 1, 2, ?)", [t1])
    DBInterface.execute(db,
        "INSERT INTO loads (id, experiment_id, load_index, start_time) VALUES (103, 1, 3, ?)", [t2])

    # Confirm they all have NULL session_id (the pre-backfill state).
    loads_pre = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, session_id FROM loads WHERE experiment_id = 1 ORDER BY id"))
    @test all(ismissing(l.session_id) || l.session_id === nothing for l in loads_pre)

    # Re-run the backfill (the sentinel was already recorded by open_db, so we
    # need to call the function directly, bypassing the sentinel guard).
    # To test the function itself, delete the sentinel and re-run migrate_schema!.
    DBInterface.execute(db,
        "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_LOAD_SESSIONS_BACKFILL])
    HimalayaUI.backfill_load_sessions!(db)

    loads_post = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, session_id FROM loads WHERE experiment_id = 1 ORDER BY id"))
    @test all(!ismissing(l.session_id) && l.session_id !== nothing for l in loads_post)
    session_ids = Int[Int(l.session_id) for l in loads_post]
    @test session_ids[1] == 1   # t0 → session 1
    @test session_ids[2] == 1   # t1 only 30 min later → still session 1
    @test session_ids[3] == 2   # t2 is 4 h after t1 → new session 2
    @test length(unique(session_ids)) == 2

    # Idempotency: re-running the backfill does NOT re-run (sentinel guards it).
    HimalayaUI.backfill_load_sessions!(db)   # sentinel present → no-op
    loads_idem = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, session_id FROM loads WHERE experiment_id = 1 ORDER BY id"))
    @test Int[Int(l.session_id) for l in loads_idem] == session_ids

    SQLite.close(db)
end
