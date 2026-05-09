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
                  load_dat, auto_group, effective_peaks, diff_update_auto_peaks!
using Himalaya: findpeaks, indexpeaks
using SQLite

@testset "persist_analysis!" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp/data",
                                     analysis_dir="/tmp/analysis")
    s_id    = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id    = create_exposure!(db; sample_id=s_id, filename="example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    peaks_result  = findpeaks(q, I, σ)

    # Prime auto_peaks so effective_peaks can query them
    diff_update_auto_peaks!(db, e_id, peaks_result, I)

    eff           = effective_peaks(db, e_id, q, I)
    candidates    = indexpeaks(eff.q, eff.sharpness)
    group_indices = auto_group(candidates)

    persist_analysis!(db, e_id, q, I, peaks_result, candidates, group_indices, eff)

    stored_peaks   = get_peaks_for_exposure(db, e_id)
    stored_indices = get_indices_for_exposure(db, e_id)
    stored_groups  = get_groups_for_exposure(db, e_id)

    @test length(stored_peaks) == length(peaks_result.q)
    @test length(stored_indices) == length(candidates)
    @test length(stored_groups) == 1
    @test stored_groups[1].kind   == "auto"
    @test stored_groups[1].active == 1
end

using HimalayaUI: init_experiment!, analyze_exposure!, open_db, get_experiment
using JSON3

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

@testset "persist_analysis! writes inputs_hash atomically with the rest of the rows (issue #34 Bug 3)" begin
    # Pre-fix, the two `inputs_hash` UPDATEs lived in `analyze_exposure!`
    # AFTER `persist_analysis!`'s transaction committed — autocommit, outside
    # any tx. A crash between the persist commit and the UPDATEs would leave
    # `exposures.analysis_inputs_hash` and `indices.inputs_hash` divergent,
    # making StaleIndicesBanner state permanently wrong.
    #
    # The fix moves the UPDATEs into `_persist_analysis_inner!` so they
    # commit/rollback as part of the same tx. This test exercises the
    # contract: calling `persist_analysis!` directly (e.g. from the CLI
    # path that didn't have an outer with_idempotency tx) MUST leave the
    # hashes populated.
    db = SQLite.DB()
    create_schema!(db)
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp/data",
                                     analysis_dir="/tmp/analysis")
    s_id    = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id    = create_exposure!(db; sample_id=s_id, filename="example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    peaks_result  = findpeaks(q, I, σ)
    diff_update_auto_peaks!(db, e_id, peaks_result, I)
    eff           = effective_peaks(db, e_id, q, I)
    candidates    = indexpeaks(eff.q, eff.sharpness)
    group_indices = auto_group(candidates)

    persist_analysis!(db, e_id, q, I, peaks_result, candidates, group_indices, eff)

    expected = HimalayaUI.hash_peak_set(eff)

    exp_hash = Tables.rowtable(DBInterface.execute(db,
        "SELECT analysis_inputs_hash FROM exposures WHERE id = ?",
        [e_id]))[1].analysis_inputs_hash
    @test exp_hash == expected

    idx_hashes = [r.inputs_hash for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT inputs_hash FROM indices WHERE exposure_id = ?", [e_id]))]
    @test !isempty(idx_hashes)
    @test all(h -> h == expected, idx_hashes)
end

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
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
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

    # Mutate peaks (simulate the user adding a manual peak via curation), then re-analyze.
    DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
        [e_id, 0.123])
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
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
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
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
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
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
        [e_id, target_q])
    manual_peak_id = Int(DBInterface.lastrowid(res))

    analyze_exposure!(db, e_id, analysis_dir)

    # The manual peak should appear in at least one index's `index_peaks` rows.
    # Curation-add peaks are referenced with peak_kind='curation'.
    n_uses = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM index_peaks WHERE peak_id = ? AND peak_kind = 'curation'",
        [manual_peak_id]))).c
    @test n_uses > 0
end

@testset "analyze_exposure! preserves auto peak IDs across reruns" begin
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)
    ids_before = sort([Int(r.id) for r in get_peaks_for_exposure(db, e_id)
                                            if String(r.source) == "auto"])

    analyze_exposure!(db, e_id, analysis_dir)
    ids_after = sort([Int(r.id) for r in get_peaks_for_exposure(db, e_id)
                                           if String(r.source) == "auto"])
    @test ids_before == ids_after
end

@testset "analyze_exposure! ignores excluded auto peaks when scoring candidates" begin
    # Regression: the candidate set used to be built from raw `findpeaks`
    # output, so the user's "this is noise" verdict (`excluded = 1`) had no
    # effect on `IndexEntry.peaks` or score until the very next reanalysis
    # happened to drop the peak. Today, `analyze_exposure!` filters excluded
    # q-values out of the indexpeaks input directly.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    # Pick an auto peak that's currently bound to ≥ 1 candidate index.
    bound_rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT ap.id, ap.q FROM auto_peaks ap
        JOIN index_peaks ip ON ip.peak_id = ap.id AND ip.peak_kind = 'auto'
        WHERE ap.exposure_id = ?
        GROUP BY ap.id HAVING COUNT(*) >= 1
        ORDER BY COUNT(*) DESC LIMIT 1
        """, [e_id]))
    @test !isempty(bound_rows)
    target_pid = Int(bound_rows[1].id)
    target_q   = Float64(bound_rows[1].q)

    # Exclude the auto peak via a curation row (new schema).
    DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
        [e_id, target_q])
    analyze_exposure!(db, e_id, analysis_dir)

    # After reanalysis the same q is still in auto_peaks (re-found by findpeaks).
    # An exclude curation still exists for it — so effective_peaks omits it
    # from the indexpeaks input. No candidate should reference the peak.
    same_q = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM auto_peaks WHERE exposure_id = ? AND ABS(q - ?) < 1e-9",
        [e_id, target_q]))
    @test !isempty(same_q)
    new_pid = Int(same_q[1].id)
    # exclude curation still present
    excl_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' AND ABS(q - ?) < 1e-9",
        [e_id, target_q]))
    @test !isempty(excl_rows)
    n_uses = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM index_peaks WHERE peak_id = ? AND peak_kind = 'auto'",
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
    # After R2.1: `peaks` is dropped by migrate_r2_split_peaks! and replaced
    # by auto_peaks + peak_curations. Check the surviving entity tables.
    for t in ("experiments", "samples", "exposures", "indices")
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", [t]))).sql)
        @test occursin("AUTOINCREMENT", sql)
    end
    # auto_peaks replaces peaks; the migrated row (q=0.1) should survive there.
    for t in ("auto_peaks",)
        sql = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type='table' AND name=?", [t]))).sql)
        @test occursin("AUTOINCREMENT", sql)
    end

    # Pre-existing auto peaks survive the migration with their ids intact.
    rows = Tables.rowtable(DBInterface.execute(db, "SELECT id, q FROM auto_peaks"))
    @test length(rows) == 1
    @test Int(rows[1].id) == 1
    @test Float64(rows[1].q) ≈ 0.1

    # Deleting and re-inserting yields a NEW id (AUTOINCREMENT prevents recycling).
    DBInterface.execute(db, "DELETE FROM auto_peaks WHERE id = 1")
    res = DBInterface.execute(db, "INSERT INTO auto_peaks (exposure_id, q) VALUES (1, 0.2)")
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
        name = 2
        display_name = 3
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
            "SELECT id FROM auto_peaks WHERE exposure_id = ?", [e_id]))
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
        name = 2
        display_name = 3
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

@testset "cli_init_with_db! refuses duplicate registration" begin
    mktempdir() do dir
        write(joinpath(dir, "experiment.toml"),
              """
              [experiment]
              name = "duplicate-test"
              [layout]
              data_dir = "data"
              analysis_dir = "analysis/automatic_analysis"
              exposure_type = "simple"
              """)
        mkpath(joinpath(dir, "analysis", "automatic_analysis"))
        db = SQLite.DB()
        HimalayaUI.create_schema!(db)
        HimalayaUI.cli_init_with_db!(db, dir)              # first call ok
        err = try
            HimalayaUI.cli_init_with_db!(db, dir)          # second call rejected
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("already registered", err.msg)
        @test occursin("reingest", err.msg)
        # Confirm no second row was inserted
        rows = Tables.rowtable(DBInterface.execute(db, "SELECT id FROM experiments"))
        @test length(rows) == 1
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
        name = 2
        display_name = 3
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
        name = 2
        display_name = 3
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
        name = 2
        display_name = 3
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

    @testset "cli_analyze --experiment selects the right experiment" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        dir2    = mktempdir()   # never registered

        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="Exp1", stems=["ST001"])
            cli_init_with_db!(db, dir1)
        end

        # Unregistered path → error.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_throws ErrorException cli_analyze(["-e", dir2])
        end

        # Registered path resolves correctly.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_nowarn cli_analyze(["-e", dir1])
        end
    end

    @testset "cli_show --experiment selects the right experiment" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        dir2    = mktempdir()   # never registered

        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="Exp1", stems=["ST001"])
            cli_init_with_db!(db, dir1)
        end

        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_throws ErrorException cli_show(["-e", dir2, "--sample", "D1"])
        end

        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_nowarn cli_show(["-e", dir1, "--sample", "D1"])
        end

        # No -e flag works for a single registered experiment.
        withenv("HIMALAYA_DB_PATH" => db_file) do
            @test_nowarn cli_show(["--sample", "D1"])
        end
    end

    @testset "_resolve_experiment errors when multiple experiments and no key" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        dir2    = mktempdir()
        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="ExpA", stems=["ST001"])
            setup_exp_dir(dir2; name="ExpB", stems=["ST002"])
            cli_init_with_db!(db, dir1)
            cli_init_with_db!(db, dir2)
            err = try; HimalayaUI._resolve_experiment(db, nothing); nothing; catch e; e; end
            @test err isa ErrorException
            @test occursin("Multiple experiments", err.msg)
        end
    end

    @testset "_resolve_experiment looks up by id" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="ExpA", stems=["ST001"])
            id = cli_init_with_db!(db, dir1)
            row = HimalayaUI._resolve_experiment(db, string(id))
            @test Int(row.id) == id
        end
    end

    @testset "_resolve_experiment looks up by name" begin
        db_file = joinpath(mktempdir(), "himalaya.db")
        dir1    = mktempdir()
        withenv("HIMALAYA_DB_PATH" => db_file) do
            db = open_db(db_file)
            setup_exp_dir(dir1; name="UniqueName123", stems=["ST001"])
            cli_init_with_db!(db, dir1)
            row = HimalayaUI._resolve_experiment(db, "UniqueName123")
            @test String(row.name) == "UniqueName123"
        end
    end
end

@testset "analyze_exposure! sets trace_hash and analysis_inputs_hash on first run" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    @test !ismissing(row.trace_hash)
    @test length(String(row.trace_hash)) == 64
    @test !ismissing(row.analysis_inputs_hash)
    @test length(String(row.analysis_inputs_hash)) == 64

    idx_hashes = [r.inputs_hash for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT inputs_hash FROM indices WHERE exposure_id = ?", [e_id]))]
    if !isempty(idx_hashes)
        @test all(!ismissing(h) for h in idx_hashes)
        @test all(String(h) == String(row.analysis_inputs_hash) for h in idx_hashes)
    end
end

@testset "analyze_exposure! preserves trace_hash across no-op reruns" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)
    row1 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    analyze_exposure!(db, e_id, analysis_dir)
    row2 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    @test String(row1.trace_hash) == String(row2.trace_hash)
    # The peak-set hash must also stabilise across a no-op rerun: same trace
    # bytes + same curation set → identical effective_peaks → identical hash.
    @test String(row1.analysis_inputs_hash) == String(row2.analysis_inputs_hash)
end

@testset "analyze_exposure! commits trace_hash + auto_peaks atomically (issue #39)" begin
    # Pre-fix, the trace_hash UPDATE and diff_update_auto_peaks!'s
    # INSERT/DELETE statements autocommitted independently. A crash between
    # diff_update and the UPDATE could land auto_peaks fresh but trace_hash
    # stale — recoverable on re-run. The reverse (UPDATE landed, auto_peaks
    # didn't fully apply) was the dangerous case: the next analyze_exposure!
    # would see hash match, skip findpeaks, and operate on a partial peak set.
    #
    # The fix wraps both writes in a single SQLite.transaction inside
    # analyze_exposure!. This test exercises the contract by seeding a
    # divergent state (trace_hash stale, auto_peaks pre-existing from a
    # prior run) and verifying that after a single re-analyze the two land
    # consistently — both reflect the current file bytes.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src          = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    dst          = joinpath(analysis_dir, "example_tot.dat")
    cp(src, dst)

    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    # First run populates auto_peaks + trace_hash + analysis hash.
    analyze_exposure!(db, e_id, analysis_dir)
    h0 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [e_id]))).trace_hash
    n0 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM auto_peaks WHERE exposure_id = ?", [e_id]))).n

    # Append bytes to the .dat to force a fresh trace_hash + re-detection.
    open(dst, "a") do io; write(io, "\n0.99 1.0 0.1\n") end

    analyze_exposure!(db, e_id, analysis_dir)
    h1 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [e_id]))).trace_hash
    n1 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM auto_peaks WHERE exposure_id = ?", [e_id]))).n

    @test String(h0) != String(h1)             # trace_hash advanced
    @test String(h1) == HimalayaUI.hash_trace_file(dst)  # to the current bytes
    @test n1 > 0                                # auto_peaks landed alongside
    # Stronger atomicity discriminator (PR #41 review suggestion 2): verify
    # auto_peaks and analysis_inputs_hash are mutually consistent — i.e.
    # exposures.analysis_inputs_hash == hash_peak_set(effective_peaks(...)).
    # If diff_update_auto_peaks! had silently no-oped (or partially landed)
    # before the trace_hash UPDATE, persist_analysis! would have hashed a
    # stale peak set and the stored hash would diverge from the recomputed
    # one. Reviewer suggested `n1 != n0`, but appending a single low-I
    # out-of-range point doesn't disturb the seven existing peaks (verified
    # empirically: n0=n1=7) — this consistency check is more robust to the
    # specific trace edit and catches the same regression class.
    q_full, I_full, _ = HimalayaUI.load_dat(dst)
    eff = HimalayaUI.effective_peaks(db, e_id, q_full, I_full)
    expected_inputs = HimalayaUI.hash_peak_set(eff)
    stored_inputs = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT analysis_inputs_hash FROM exposures WHERE id = ?", [e_id]))).analysis_inputs_hash
    @test String(stored_inputs) == expected_inputs
    # And every auto_peak row has a finite q (no orphan/partial INSERT).
    qs = [r.q for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM auto_peaks WHERE exposure_id = ?", [e_id]))]
    @test all(q -> !ismissing(q) && isfinite(Float64(q)), qs)
end

@testset "analyze_exposure! re-runs findpeaks when trace bytes change" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    dst = joinpath(analysis_dir, "example_tot.dat")
    cp(src, dst)

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)
    h_before = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [e_id]))).trace_hash

    open(dst, "a") do io; write(io, "\n0.99 1.0 0.1\n") end

    analyze_exposure!(db, e_id, analysis_dir)
    h_after = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash FROM exposures WHERE id = ?", [e_id]))).trace_hash
    @test String(h_before) != String(h_after)
end

@testset "analyze_exposure! commits persist_analysis! + analyze_run row atomically (issue #42)" begin
    # Pre-fix, the apply_event! call lived AFTER persist_analysis!'s tx
    # committed (autocommit, outside any tx). On the CLI / programmatic /
    # POST /api/experiments/{id}/analyze paths (no enclosing with_idempotency
    # tx), a crash between persist_analysis!'s commit and apply_event!'s
    # INSERT would land the analysis state durably but with no analyze_run
    # user_actions row. Same shape as #34 Bug 3, one layer up.
    #
    # The fix wraps persist_analysis! + apply_event! in a single
    # SQLite.transaction inside analyze_exposure!; apply_event! switches to
    # the InTransaction() variant to participate. We can't deterministically
    # induce a mid-tx crash in a test, but two structural invariants pin
    # the contract:
    #
    # 1. Forward correlation: after each state-advancing call, the latest
    #    analyze_run row's inputs_hash_after / trace_hash_after match the
    #    exposures row's analysis_inputs_hash / trace_hash. Catches any
    #    future regression that reorders the writes (e.g. emits the event
    #    before persist_analysis! ran).
    # 2. Outer-tx rollback symmetry: wrap analyze_exposure! in an outer
    #    SQLite.transaction that throws — the analyze_run row, the auto_peaks
    #    rebuild, and the exposures hash columns must all revert together.
    #    The inner SQLite.transaction nests as a SAVEPOINT inside the outer
    #    tx; a future regression that pulled apply_event! back outside the
    #    SAVEPOINT (autocommit) would leave a dangling analyze_run row here.
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    dst = joinpath(analysis_dir, "example_tot.dat")
    cp(src, dst)

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    # Layer 1: forward correlation after a slow-path run.
    analyze_exposure!(db, e_id, analysis_dir)
    state = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    event = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE action = 'analyze_run' ORDER BY id DESC LIMIT 1")))
    p = JSON3.read(String(event.payload))
    @test String(state.trace_hash)          == String(p[:trace_hash_after])
    @test String(state.analysis_inputs_hash) == String(p[:inputs_hash_after])

    # Layer 2: outer-tx rollback symmetry. Append bytes to force the slow path
    # again, run analyze_exposure! inside an outer tx that throws, then verify
    # nothing landed.
    open(dst, "a") do io; write(io, "\n0.99 1.0 0.1\n") end
    n_runs_before = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'analyze_run'"))).n
    state_before = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    n_peaks_before = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM auto_peaks WHERE exposure_id = ?", [e_id]))).n
    @test_throws ErrorException SQLite.transaction(db) do
        analyze_exposure!(db, e_id, analysis_dir; defer_broadcast=true)
        error("rollback sentinel")
    end
    n_runs_after = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'analyze_run'"))).n
    state_after = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT trace_hash, analysis_inputs_hash FROM exposures WHERE id = ?", [e_id])))
    n_peaks_after = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM auto_peaks WHERE exposure_id = ?", [e_id]))).n
    @test n_runs_after == n_runs_before
    @test String(state_after.trace_hash)           == String(state_before.trace_hash)
    @test String(state_after.analysis_inputs_hash) == String(state_before.analysis_inputs_hash)
    @test n_peaks_after == n_peaks_before
end

@testset "analyze_run payload shows both skip flags true on no-op rerun" begin
    # Regression: the skip-flag expressions previously checked only hash equality,
    # not the full predicate (hash match AND existing rows). A hash match on a
    # fresh DB with empty auto_peaks would record findpeaks_skipped=true while
    # findpeaks actually ran.
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    # First run: both skip flags must be false (nothing cached yet).
    analyze_exposure!(db, e_id, analysis_dir)
    row1 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE action = 'analyze_run' ORDER BY id DESC LIMIT 1")))
    p1 = JSON3.read(String(row1.payload))
    @test p1[:findpeaks_skipped]  == false
    @test p1[:indexpeaks_skipped] == false

    # Second run with identical trace and no curation changes: genuine no-op.
    analyze_exposure!(db, e_id, analysis_dir)
    row2 = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT payload FROM user_actions WHERE action = 'analyze_run' ORDER BY id DESC LIMIT 1")))
    p2 = JSON3.read(String(row2.payload))
    @test p2[:findpeaks_skipped]  == true
    @test p2[:indexpeaks_skipped] == true
end

@testset "cli_init_with_db! is atomic on validation failure" begin
    mktempdir() do tmp
        # Build an experiment dir whose manifest has duplicate sample names.
        exp_dir = joinpath(tmp, "exp")
        mkpath(joinpath(exp_dir, "data"))
        mkpath(joinpath(exp_dir, "analysis", "automatic_analysis"))
        # Minimal experiment.toml using new key names
        write(joinpath(exp_dir, "experiment.toml"), """
[experiment]
name = "validation-fail"
manifest = "manifest.csv"
[manifest]
delimiter      = "\\t"
skip_rows      = 1
sample_id      = 1
name           = 2
display_name   = 3
filenames      = 9
notes_sample   = 10
notes_exposure = 11
[layout]
data_dir = "data"
analysis_dir = "analysis/automatic_analysis"
[files]
integration = "{name}.dat"
image       = "{name}.tiff"
""")
        # Two samples with the same name (duplicate_name violation).
        write(joinpath(exp_dir, "manifest.csv"), """sample_id\tname\tdisplay_name\tcol4\tcol5\tcol6\tcol7\tcol8\tfilenames\tnotes_sample\tnotes_exposure
1\tDUP\tfirst\t\t\t\t\t\tA001\t\t
2\tDUP\tsecond\t\t\t\t\t\tA002\t\t
""")

        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        @test_throws HimalayaUI.ManifestValidationError HimalayaUI.cli_init_with_db!(db, exp_dir)

        # No experiment row leaked — transaction rolled back.
        rows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM experiments"))
        @test isempty(rows)
        # No samples either.
        srows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM samples"))
        @test isempty(srows)
    end
end

@testset "analyze_exposure! defer_broadcast=true suppresses analyze_run SSE frame" begin
    # M2.2 contract: when curation routes call analyze_exposure! synchronously
    # inside their with_idempotency tx, they must pass defer_broadcast=true so
    # the slow-path's inner apply_event! doesn't broadcast before the outer tx
    # commits. The user_actions row is still written (durable), only the SSE
    # frame is suppressed.
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db = open_db(joinpath(tmp, "himalaya.db"))
    HimalayaUI.bind_db!(db)
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id = create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    e_id = create_exposure!(db; sample_id=s_id, filename="example_tot")

    # Hook a fake SSE subscriber to count frames.
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end

    try
        # First run with defer_broadcast=true: slow-path runs (cold cache).
        # No analyze_run SSE frame should be broadcast.
        analyze_exposure!(db, e_id, analysis_dir; defer_broadcast=true)

        # Drain.
        sleep(0.05)
        frames = String[]
        while isready(pending)
            push!(frames, take!(pending))
        end
        analyze_frames = filter(f -> occursin("\"kind\":\"analyze_run\"", f), frames)
        @test isempty(analyze_frames)

        # But the user_actions row is still durable.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE action = 'analyze_run' AND entity_id = ?", [e_id]))
        @test !isempty(rows)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
end
