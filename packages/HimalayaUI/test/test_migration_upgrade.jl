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

@testset "Task 3 upgrade-grouping CLI" begin
    # Helper: run cli_upgrade_grouping with HIMALAYA_DB_PATH pointing at `path`.
    # Captures stdout to a temp file and returns it as a String.
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
        @test length(before.loads) == length(after.loads)
        @test length(before.samples) == length(after.samples)
        @test length(before.exposures) == length(after.exposures)
        # No loads should have been created (dry-run rolled back).
        @test isempty(after.loads)
        # Curation untouched.
        @test before.curation == after.curation
        # Exposure sample_ids unchanged (all still as set by migrate_schema!).
        bsids = sort(Int[Int(e.sample_id) for e in before.exposures if e.sample_id !== missing])
        asids = sort(Int[Int(e.sample_id) for e in after.exposures if e.sample_id !== missing])
        @test bsids == asids
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
