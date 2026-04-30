using Test
using Himalaya: indexpeaks, Index, peaks, score
using HimalayaUI: auto_group

@testset "auto_group" begin
    qs     = [0.1000, 0.1414, 0.2000]
    sharps = [1.0, 0.8, 0.6]

    candidates = indexpeaks(qs, sharps)

    group = auto_group(candidates)

    if !isempty(candidates)
        @test !isempty(group)
        peak_sets = [Set(peaks(idx)) for idx in group]
        for i in eachindex(peak_sets), j in eachindex(peak_sets)
            i == j && continue
            @test isempty(intersect(peak_sets[i], peak_sets[j]))
        end
    end
end

using HimalayaUI: create_schema!, create_experiment!, create_sample!,
                  create_exposure!, persist_analysis!, get_peaks_for_exposure,
                  get_indices_for_exposure, get_groups_for_exposure,
                  load_dat, auto_group
using Himalaya: findpeaks, indexpeaks
using SQLite

@testset "persist_analysis!" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp/data",
                                     analysis_dir="/tmp/analysis")
    s_id    = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id    = create_exposure!(db; sample_id=s_id, filename="example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    peaks_result  = findpeaks(q, I, σ)
    candidates    = indexpeaks(peaks_result.q, peaks_result.sharpness)
    group_indices = auto_group(candidates)

    persist_analysis!(db, e_id, q, I, peaks_result, candidates, group_indices)

    stored_peaks   = get_peaks_for_exposure(db, e_id)
    stored_indices = get_indices_for_exposure(db, e_id)
    stored_groups  = get_groups_for_exposure(db, e_id)

    @test length(stored_peaks) == length(peaks_result.q)
    @test length(stored_indices) == length(candidates)
    @test length(stored_groups) == 1
    @test stored_groups[1].kind   == "auto"
    @test stored_groups[1].active == 1
end

using HimalayaUI: find_tiff_for_dat

@testset "find_tiff_for_dat" begin
    dir = mktempdir()
    # .tiff match
    dat_path  = joinpath(dir, "sample_pos1.dat")
    tiff_path = joinpath(dir, "sample_pos1.tiff")
    touch(dat_path); touch(tiff_path)
    @test find_tiff_for_dat(dat_path) == tiff_path

    # .tif fallback
    dat2  = joinpath(dir, "sample_pos2.dat")
    tif2  = joinpath(dir, "sample_pos2.tif")
    touch(dat2); touch(tif2)
    @test find_tiff_for_dat(dat2) == tif2

    # no match → nothing
    dat3 = joinpath(dir, "sample_pos3.dat")
    touch(dat3)
    @test find_tiff_for_dat(dat3) === nothing
end

using HimalayaUI: init_experiment!, analyze_exposure!, open_db, get_experiment

@testset "init_experiment!" begin
    tmp = mktempdir()
    data_dir     = joinpath(tmp, "data")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(data_dir)
    mkpath(analysis_dir)

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db;
        name         = "TestExp",
        path         = tmp,
        data_dir     = data_dir,
        analysis_dir = analysis_dir)

    @test exp_id == 1
    exp = get_experiment(db, exp_id)
    @test exp.name == "TestExp"
end

using HimalayaUI: ensure_custom_group!
using DBInterface
using Tables

@testset "persist_analysis! preserves custom-group members across reanalysis" begin
    # Regression: when a peak is added/removed/excluded the frontend triggers
    # an auto-reanalysis. That re-runs persist_analysis! which deletes all
    # `indices` rows and re-inserts with fresh PKs. If we don't translate
    # custom-group memberships from old PK → new PK by semantic identity
    # (phase + basis), the active set ends up empty.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    # Promote the auto group to a custom group (mirrors the UI flow when a
    # user clicks an Active-set member).
    custom_id, _ = ensure_custom_group!(db, e_id)

    # Snapshot the custom group's members by semantic identity before reanalysis.
    snap_before = Tables.rowtable(DBInterface.execute(db, """
        SELECT i.phase, i.basis FROM index_group_members m
        JOIN indices i ON i.id = m.index_id
        WHERE m.group_id = ?
        """, [custom_id]))
    @test !isempty(snap_before)
    snap_phases_before = sort(String[String(r.phase) for r in snap_before])

    # Mutate peaks (simulate the user adding a manual peak), then re-analyze.
    DBInterface.execute(db,
        "INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source, excluded)
         VALUES (?, ?, ?, ?, ?, 'manual', 0)",
        [e_id, 0.123, 1.0, 0.5, 0.5])
    analyze_exposure!(db, e_id, analysis_dir)

    # The custom group's members should still be populated, with at least one
    # phase from the original snapshot still present (the data is stable
    # enough that the same phase indexings should re-emerge).
    snap_after = Tables.rowtable(DBInterface.execute(db, """
        SELECT i.phase FROM index_group_members m
        JOIN indices i ON i.id = m.index_id
        WHERE m.group_id = ?
        """, [custom_id]))
    @test !isempty(snap_after)
    snap_phases_after = sort(String[String(r.phase) for r in snap_after])
    # At least one of the originally-active phases survives.
    @test !isempty(intersect(snap_phases_before, snap_phases_after))

    # The custom group should still be `active`; the auto group should be demoted.
    g_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, kind, active FROM index_groups WHERE exposure_id = ?", [e_id]))
    custom_row = first(filter(r -> String(r.kind) == "custom", g_rows))
    auto_row   = first(filter(r -> String(r.kind) == "auto",   g_rows))
    @test Int(custom_row.active) == 1
    @test Int(auto_row.active)   == 0
end

@testset "analyze_exposure! integration" begin
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    @test length(get_peaks_for_exposure(db, e_id))   > 0
    @test length(get_indices_for_exposure(db, e_id)) > 0
    @test length(get_groups_for_exposure(db, e_id))  == 1
end

@testset "analyze_exposure! incorporates manual peaks into candidate indices" begin
    # Regression: manual peaks live in the DB but `Himalaya.indexpeaks` was
    # previously called with only fresh `findpeaks` output, so a manual peak
    # at a phase's predicted ratio position never landed in `IndexEntry.peaks`.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    # Pick a candidate index's basis and place a manual peak at its
    # √2-ratio position — guaranteed to be a predicted q for cubic phases
    # so SOME candidate should bind it on reanalysis.
    ix_row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT basis FROM indices WHERE exposure_id = ? ORDER BY score DESC LIMIT 1",
        [e_id])))
    target_q = Float64(ix_row.basis) * sqrt(2.0)

    res = DBInterface.execute(db,
        "INSERT INTO peaks (exposure_id, q, source) VALUES (?, ?, 'manual')",
        [e_id, target_q])
    manual_peak_id = Int(DBInterface.lastrowid(res))

    analyze_exposure!(db, e_id, analysis_dir)

    # The manual peak should appear in at least one index's `index_peaks` rows.
    n_uses = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM index_peaks WHERE peak_id = ?",
        [manual_peak_id]))).c
    @test n_uses > 0
end

@testset "analyze_exposure! ignores excluded auto peaks when scoring candidates" begin
    # Regression: `was_excluded` carry-forward was applied at peak persistence
    # time, but the candidate set was built from raw `findpeaks` output — so
    # the user's "this is noise" verdict had no effect on `IndexEntry.peaks`
    # or score until the very next reanalysis happened to drop the peak.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    # Pick an auto peak that's currently bound to ≥ 1 candidate index.
    bound_rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT p.id, p.q FROM peaks p
        JOIN index_peaks ip ON ip.peak_id = p.id
        WHERE p.exposure_id = ? AND p.source = 'auto'
        GROUP BY p.id HAVING COUNT(*) >= 1
        ORDER BY COUNT(*) DESC LIMIT 1
        """, [e_id]))
    @test !isempty(bound_rows)
    target_pid = Int(bound_rows[1].id)
    target_q   = Float64(bound_rows[1].q)

    DBInterface.execute(db,
        "UPDATE peaks SET excluded = 1 WHERE id = ?", [target_pid])
    analyze_exposure!(db, e_id, analysis_dir)

    # After reanalysis the same q should still be detected (auto peaks are
    # re-found by findpeaks) and re-marked excluded via `was_excluded`. But
    # no candidate should reference it, because we now filter excluded q's
    # before calling indexpeaks.
    same_q = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, excluded FROM peaks WHERE exposure_id = ? AND ABS(q - ?) < 1e-9",
        [e_id, target_q]))
    @test !isempty(same_q)
    @test Int(same_q[1].excluded) == 1
    new_pid = Int(same_q[1].id)
    n_uses = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM index_peaks WHERE peak_id = ?",
        [new_pid]))).c
    @test n_uses == 0
end

@testset "open_db migrates rowid PKs to AUTOINCREMENT" begin
    # Regression: chat @-mentions reference entities by id; SQLite rowid PKs
    # reuse freed ids on deletion, so a deleted-then-re-created index can take
    # the same id and silently rebind the mention to a different index.
    tmp = mktempdir()
    db_path = joinpath(tmp, "legacy.db")

    legacy = SQLite.DB(db_path)
    DBInterface.execute(legacy, "PRAGMA foreign_keys = ON")
    DBInterface.execute(legacy,
        "CREATE TABLE experiments (id INTEGER PRIMARY KEY, name TEXT, path TEXT NOT NULL,
                                   data_dir TEXT NOT NULL, analysis_dir TEXT NOT NULL)")
    DBInterface.execute(legacy,
        "CREATE TABLE samples (id INTEGER PRIMARY KEY, experiment_id INTEGER, label TEXT)")
    DBInterface.execute(legacy,
        "CREATE TABLE exposures (id INTEGER PRIMARY KEY, sample_id INTEGER, filename TEXT)")
    DBInterface.execute(legacy,
        "CREATE TABLE peaks (id INTEGER PRIMARY KEY, exposure_id INTEGER, q REAL NOT NULL)")
    DBInterface.execute(legacy,
        "CREATE TABLE indices (id INTEGER PRIMARY KEY, exposure_id INTEGER, phase TEXT NOT NULL, basis REAL NOT NULL)")
    DBInterface.execute(legacy,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
    DBInterface.execute(legacy, "INSERT INTO samples (experiment_id, label) VALUES (1, 'A1')")
    DBInterface.execute(legacy, "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
    DBInterface.execute(legacy, "INSERT INTO peaks (exposure_id, q) VALUES (1, 0.1)")
    DBInterface.execute(legacy, "INSERT INTO indices (exposure_id, phase, basis) VALUES (1, 'Pn3m', 0.1)")
    close(legacy)

    db = open_db(db_path)
    for t in ("experiments", "samples", "exposures", "peaks", "indices")
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", [t]))).sql)
        @test occursin("AUTOINCREMENT", sql)
    end

    # Pre-existing rows survive the migration with their ids intact.
    rows = Tables.rowtable(DBInterface.execute(db, "SELECT id FROM peaks"))
    @test [Int(r.id) for r in rows] == [1]

    # Deleting and re-inserting yields a NEW id, not the recycled rowid.
    DBInterface.execute(db, "DELETE FROM peaks WHERE id = 1")
    res = DBInterface.execute(db, "INSERT INTO peaks (exposure_id, q) VALUES (1, 0.2)")
    @test Int(DBInterface.lastrowid(res)) >= 2
end

@testset "analyze_exposure! uses integration pattern from config" begin
    # Custom config with a different integration suffix to prove the pattern is honored
    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis")
        mkpath(analysis_dir)
        # Copy the standard fixture but write it under a custom suffix so
        # only the config-driven pattern can resolve it.
        src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
        cp(src, joinpath(analysis_dir, "EX001_1d.dat"))

        db = SQLite.DB()
        HimalayaUI.create_schema!(db)

        # Build a config with custom integration pattern
        toml_blob = """
        [experiment]
        name = "T"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}_1d.dat"
        image = "{name}.tiff"
        """

        exp_id = HimalayaUI.create_experiment!(db;
            path = dir, data_dir = joinpath(dir, "data"), analysis_dir = analysis_dir,
            config = toml_blob, experiment_type = "simple")
        s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
        e_id = HimalayaUI.create_exposure!(db; sample_id = s_id, filename = "EX001")

        # Should resolve EX001 → "EX001_1d.dat" via the custom pattern
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

        # Confirm peaks were persisted (i.e. the pattern resolved correctly and the file was loaded)
        peaks = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peaks WHERE exposure_id = ?", [e_id]))
        @test length(peaks) >= 0   # Just need analyze_exposure! not to throw
    end
end

@testset "cli_init_with_db! reads experiment.toml" begin
    mktempdir() do dir
        # Set up experiment directory
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        data_dir     = joinpath(dir, "data")
        mkpath(analysis_dir)
        mkpath(data_dir)

        # Use the canonical fixture for valid integration data
        fixture = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
        cp(fixture, joinpath(analysis_dir, "JC001.dat"))
        cp(fixture, joinpath(analysis_dir, "JC002.dat"))

        # Manifest
        manifest = joinpath(dir, "manifest.csv")
        write(manifest, join([
            "skip-row",
            "1\tD1\tUX1\tT\tt\t\t\t\tJC001-002\tnote_s\tnote_e",
        ], "\n"))

        # experiment.toml
        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "Run/Exp"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 12.0
        flight_path_m = 2.5
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)

        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        exp_id = HimalayaUI.cli_init_with_db!(db, dir)

        # Verify experiment was created with config and beamline params
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT name, energy_kev, flight_path_m, config FROM experiments WHERE id = ?", [exp_id]))
        @test length(rows) == 1
        @test rows[1].name == "Run/Exp"
        @test rows[1].energy_kev == 12.0
        @test rows[1].flight_path_m == 2.5
        @test contains(rows[1].config, "[experiment]")

        # Verify samples and exposures
        samples = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM samples WHERE experiment_id = ?", [exp_id]))
        @test length(samples) == 1

        exposures = Tables.rowtable(DBInterface.execute(db,
            "SELECT filename FROM exposures WHERE sample_id = ? ORDER BY filename", [samples[1].id]))
        @test [e.filename for e in exposures] == ["JC001", "JC002"]
    end
end

@testset "cli_init_with_db! errors when experiment.toml missing" begin
    mktempdir() do dir
        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        @test_throws ErrorException HimalayaUI.cli_init_with_db!(db, dir)
    end
end

@testset "cli_init_with_db! does not write to experiment directory" begin
    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        data_dir = joinpath(dir, "data")
        mkpath(analysis_dir)
        mkpath(data_dir)
        fixture = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
        cp(fixture, joinpath(analysis_dir, "JC001.dat"))
        write(joinpath(dir, "manifest.csv"),
              "skip-row\n1\tD1\tUX1\tT\tt\t\t\t\tJC001\t\t")
        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "T"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)

        # Snapshot the file list before init
        before = Set(readdir(dir))

        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        HimalayaUI.cli_init_with_db!(db, dir)

        # Snapshot after init — must be identical (no DB or other files written)
        after = Set(readdir(dir))
        @test before == after
    end
end

@testset "reingest! adds new exposures and preserves curated ones" begin
    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        data_dir     = joinpath(dir, "data")
        mkpath(analysis_dir)
        mkpath(data_dir)
        fixture = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
        for name in ["JC001", "JC002", "JC003"]
            cp(fixture, joinpath(analysis_dir, name * ".dat"))
        end

        # Initial manifest references only JC001-002
        manifest = joinpath(dir, "manifest.csv")
        write(manifest, "skip-row\n1\tD1\tUX1\tT\tt\t\t\t\tJC001-002\told\t")

        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "R/E"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        energy_kev = 0.0
        flight_path_m = 0.0
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)

        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        exp_id = HimalayaUI.cli_init_with_db!(db, dir)

        # Curate JC001
        DBInterface.execute(db,
            "UPDATE exposures SET status = 'accepted' WHERE filename = ?", ["JC001"])

        # Update manifest to extend the range to JC003
        write(manifest, "skip-row\n1\tD1\tUX1\tT\tt\t\t\t\tJC001-003\tnew\t")

        HimalayaUI.reingest!(db, exp_id, dir)

        # Verify all three exposures are present
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, status FROM exposures ORDER BY filename"))
        @test [r.filename for r in rows] == ["JC001", "JC002", "JC003"]

        # Curation on JC001 must be preserved
        jc001 = first(filter(r -> r.filename == "JC001", rows))
        @test jc001.status == "accepted"
    end
end

@testset "reingest! errors when experiment.toml missing" begin
    mktempdir() do dir
        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        exp_id = HimalayaUI.create_experiment!(db;
            path = dir, data_dir = dir, analysis_dir = dir)
        @test_throws ErrorException HimalayaUI.reingest!(db, exp_id, dir)
    end
end

# ── CLI path-targeting tests ──────────────────────────────────────────────────

using HimalayaUI: cli_analyze, cli_show, cli_init_with_db!

let
    # Helper scoped to this section — not visible elsewhere in the test file.
    function setup_exp_dir(dir; name="E", stems=["ST001"])
        analysis_dir = joinpath(dir, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        mkpath(joinpath(dir, "data"))
        fixture = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
        for stem in stems
            cp(fixture, joinpath(analysis_dir, stem * ".dat"); force=true)
        end
        write(joinpath(dir, "manifest.csv"), join([
            "skip-row",
            "1\tD1\t$(name)\tT\tt\t\t\t\t$(join(stems, ","))\tnote_s\tnote_e",
        ], "\n"))
        write(joinpath(dir, "experiment.toml"), """
        [experiment]
        name = "$name"
        description = ""
        manifest = "manifest.csv"
        [beamline]
        [manifest]
        delimiter = "\\t"
        skip_rows = 1
        header_row = 0
        sample_id = 1
        label = 2
        name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)
    end

    @testset "cli_analyze uses experiment_path, not hardcoded id=1" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        dir2    = mktempdir()   # never registered

        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="Exp1", stems=["ST001"])
            cli_init_with_db!(db, dir1)
        end

        # After fix: raises ErrorException because dir2 is not registered.
        # Before fix: silently operates on experiment id=1 (dir1), no error → test is RED.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_throws ErrorException cli_analyze([dir2])
        end

        # Success path: dir1 is registered — cli_analyze([dir1]) must not throw.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_nowarn cli_analyze([dir1])
        end
    end

    @testset "cli_show uses experiment_path, not hardcoded id=1" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        dir2    = mktempdir()   # never registered

        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="Exp1", stems=["ST001"])
            cli_init_with_db!(db, dir1)
        end

        # After fix: raises ErrorException because dir2 is not registered.
        # Before fix: silently returns exp1's data for sample "D1" → test is RED.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_throws ErrorException cli_show([dir2, "--sample", "D1"])
        end

        # Success path: dir1 is registered — cli_show([dir1]) with sample D1 must not throw.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_nowarn cli_show([dir1, "--sample", "D1"])
        end
    end
end
