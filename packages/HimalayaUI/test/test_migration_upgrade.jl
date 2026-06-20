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

The fixture embeds the three shapes later tasks must exercise:

  (a) `(experiment_id, filename)` COLLISION sharing one `image_path`: two
      exposures with the same bare-stem `filename` but distinct `sample_id`
      (legal under the legacy `(sample_id, filename)` key). One carries a
      `trace_hash` + an `auto` index_groups/auto_peaks curation chain; the other
      is empty. Both point at the SAME `image_path` (the real-data dedup case).
      NOTE: in the legacy shape this collision is legal; the swap to the new
      key is Task 1's job, so `build_legacy_db` deliberately leaves the two rows
      un-deduped and the migration is NOT run here.

  (b) RESHOOT: one old sample whose two stems sit at the same slot position but
      are separated by a large timestamp gap — so a re-derivation lands them in
      TWO different loads.

  (c) MANY-OLD -> ONE-CELL: two old samples whose single stems share both
      timestamp-load AND horizontal position — so a re-derivation lands them in
      ONE (load, slot) cell.

Returns the path plus an `info` NamedTuple of the seeded ids/stems for tests.
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
